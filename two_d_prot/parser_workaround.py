#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy

import Bio.PDB
from Bio.PDB.PDBExceptions import PDBConstructionException

class MMCIFParser(Bio.PDB.MMCIFParser):

    def _build_structure(self, structure_id):

        # two special chars as placeholders in the mmCIF format
        # for item values that cannot be explicitly assigned
        # see: pdbx/mmcif syntax web page
        _unassigned = {".", "?"}

        mmcif_dict = self._mmcif_dict

        atom_serial_list = mmcif_dict["_atom_site.id"]
        atom_id_list = mmcif_dict["_atom_site.label_atom_id"]
        residue_id_list = mmcif_dict["_atom_site.label_comp_id"]
        try:
            element_list = mmcif_dict["_atom_site.type_symbol"]
        except KeyError:
            element_list = None
        ###chain_id_list = mmcif_dict["_atom_site.auth_asym_id"]
        chain_id_list = mmcif_dict["_atom_site.label_asym_id"]
        x_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_x"]]
        y_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_y"]]
        z_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_z"]]
        alt_list = mmcif_dict["_atom_site.label_alt_id"]

        # Some values are missing in "ubertemplate" CIF, so provide some
        # defaults just to make parser happy
        try:
            icode_list = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]
        except KeyError:
            print("WARN: no _atom_site.pdbx_PDB_ins_code in mmcif_dict for %s" % structure_id)
            icode_list = ['?' for i in alt_list]
        try:
            b_factor_list = mmcif_dict["_atom_site.B_iso_or_equiv"]
        except KeyError:
            print("WARN: no _atom_site.B_iso_or_equiv in mmcif_dict for %s" % structure_id)
            b_factor_list = ['0.0' for i in alt_list]
        try:
            occupancy_list = mmcif_dict["_atom_site.occupancy"]
        except KeyError:
            print("WARN: no _atom_site.occupancy in mmcif_dict for %s" % structure_id)
            occupancy_list = ['1.00' for i in alt_list]

        fieldname_list = mmcif_dict["_atom_site.group_PDB"]
        try:
            serial_list = [int(n) for n in mmcif_dict["_atom_site.pdbx_PDB_model_num"]]
        except KeyError:
            # No model number column
            serial_list = None
        except ValueError:
            # Invalid model number (malformed file)
            raise PDBConstructionException("Invalid model number") from None
        try:
            aniso_u11 = mmcif_dict["_atom_site_anisotrop.U[1][1]"]
            aniso_u12 = mmcif_dict["_atom_site_anisotrop.U[1][2]"]
            aniso_u13 = mmcif_dict["_atom_site_anisotrop.U[1][3]"]
            aniso_u22 = mmcif_dict["_atom_site_anisotrop.U[2][2]"]
            aniso_u23 = mmcif_dict["_atom_site_anisotrop.U[2][3]"]
            aniso_u33 = mmcif_dict["_atom_site_anisotrop.U[3][3]"]
            aniso_flag = 1
        except KeyError:
            # no anisotropic B factors
            aniso_flag = 0
        #### if auth_seq_id is present, we use this.
        #### Otherwise label_seq_id is used.
        ###if "_atom_site.auth_seq_id" in mmcif_dict:
        ###    seq_id_list = mmcif_dict["_atom_site.auth_seq_id"]
        ###else:
        ###    seq_id_list = mmcif_dict["_atom_site.label_seq_id"]
        seq_id_list = mmcif_dict["_atom_site.label_seq_id"]
        # Now loop over atoms and build the structure
        current_chain_id = None
        current_residue_id = None
        current_resname = None
        structure_builder = self._structure_builder
        structure_builder.init_structure(structure_id)
        structure_builder.init_seg(" ")
        # Historically, Biopython PDB parser uses model_id to mean array index
        # so serial_id means the Model ID specified in the file
        current_model_id = -1
        current_serial_id = -1
        for i in range(0, len(atom_id_list)):

            # set the line_counter for 'ATOM' lines only and not
            # as a global line counter found in the PDBParser()
            structure_builder.set_line_counter(i)

            # Try coercing serial to int, for compatibility with PDBParser
            # But do not quit if it fails. mmCIF format specs allow strings.
            try:
                serial = int(atom_serial_list[i])
            except ValueError:
                serial = atom_serial_list[i]
                warnings.warn(
                    "PDBConstructionWarning: "
                    "Some atom serial numbers are not numerical",
                    PDBConstructionWarning,
                )

            x = x_list[i]
            y = y_list[i]
            z = z_list[i]
            resname = residue_id_list[i]
            chainid = chain_id_list[i]
            altloc = alt_list[i]
            if altloc in _unassigned:
                altloc = " "
            try:
                int_resseq = int(seq_id_list[i])
            except ValueError:
                ###raise PDBConstructionException("Invalid seq id")
                continue
            icode = icode_list[i]
            if icode in _unassigned:
                icode = " "
            name = atom_id_list[i]
            # occupancy & B factor
            try:
                tempfactor = float(b_factor_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing B factor") from None
            try:
                occupancy = float(occupancy_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing occupancy") from None
            fieldname = fieldname_list[i]
            if fieldname == "HETATM":
                if resname == "HOH" or resname == "WAT":
                    hetatm_flag = "W"
                else:
                    hetatm_flag = "H"
            else:
                hetatm_flag = " "

            resseq = (hetatm_flag, int_resseq, icode)

            if serial_list is not None:
                # model column exists; use it
                serial_id = serial_list[i]
                if current_serial_id != serial_id:
                    # if serial changes, update it and start new model
                    current_serial_id = serial_id
                    current_model_id += 1
                    structure_builder.init_model(current_model_id, current_serial_id)
                    current_chain_id = None
                    current_residue_id = None
                    current_resname = None
            else:
                # no explicit model column; initialize single model
                structure_builder.init_model(current_model_id)

            if current_chain_id != chainid:
                current_chain_id = chainid
                structure_builder.init_chain(current_chain_id)
                current_residue_id = None
                current_resname = None

            if current_residue_id != resseq or current_resname != resname:
                current_residue_id = resseq
                current_resname = resname
                structure_builder.init_residue(resname, hetatm_flag, int_resseq, icode)

            coord = numpy.array((x, y, z), "f")
            element = element_list[i].upper() if element_list else None
            structure_builder.init_atom(
                name,
                coord,
                tempfactor,
                occupancy,
                altloc,
                name,
                serial_number=serial,
                element=element,
            )
            if aniso_flag == 1 and i < len(aniso_u11):
                u = (
                    aniso_u11[i],
                    aniso_u12[i],
                    aniso_u13[i],
                    aniso_u22[i],
                    aniso_u23[i],
                    aniso_u33[i],
                )
                mapped_anisou = [float(_) for _ in u]
                anisou_array = numpy.array(mapped_anisou, "f")
                structure_builder.set_anisou(anisou_array)
        # Now try to set the cell
        try:
            a = float(mmcif_dict["_cell.length_a"][0])
            b = float(mmcif_dict["_cell.length_b"][0])
            c = float(mmcif_dict["_cell.length_c"][0])
            alpha = float(mmcif_dict["_cell.angle_alpha"][0])
            beta = float(mmcif_dict["_cell.angle_beta"][0])
            gamma = float(mmcif_dict["_cell.angle_gamma"][0])
            cell = numpy.array((a, b, c, alpha, beta, gamma), "f")
            spacegroup = mmcif_dict["_symmetry.space_group_name_H-M"][0]
            spacegroup = spacegroup[1:-1]  # get rid of quotes!!
            if spacegroup is None:
                raise Exception
            structure_builder.set_symmetry(spacegroup, cell)
        except Exception:
            pass  # no cell found, so just ignore
