#!/usr/bin/env python3
# go through all channels in given directory and find out group which goes through the given polygon
# defined by SSEs positions
# e.g. for cytochrpmes channel3 tunnels go through polygon defined by SSEs F and G

import json
import utils_read_cif
import numpy as np

MODEL = ["1og2", "A", "00"]
#MODEL_1akd = ["1akd", "A", "00"]

def sses_model_mapping_1og2(model_annotation_file):
    # maping based on  https://sestra.ncbr.muni.cz/job/cyp-1akd
    mapping = {}
    print("model_annotation_file", model_annotation_file)
    with open(model_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = f"{MODEL[0]}{MODEL[1]}{MODEL[2]}"
        if not pcd in data.keys():
            pcd = f"{MODEL[0]}{MODEL[1]}"
        for sse in data[pcd]["secondary_structure_elements"]:
            if sse["start"] == 23 and sse["end"] == 26:
                mapping["A'"] = sse["label"]
            if sse["start"] == 31 and sse["end"] == 42:
                mapping["A"] = sse["label"]
            if sse["start"] == 45 and sse["end"] == 50:
                mapping["1_1"] = sse["label"]
            if sse["start"] == 61 and sse["end"] == 71:
                mapping["B"] = sse["label"]
            if sse["start"] == 82 and sse["end"] == 88:
                mapping["B'"] = sse["label"]
            if sse["start"] == 99 and sse["end"] == 112:
                mapping["C"] = sse["label"]
            if sse["start"] == 122 and sse["end"] == 140:
                mapping["D"] = sse["label"]
            if sse["start"] == 148 and sse["end"] == 164:
                mapping["E"] = sse["label"]
            if sse["start"] == 173 and sse["end"] == 190: #H16
                mapping["F"] = sse["label"]
            if sse["start"] == 192 and sse["end"] == 200: #H20
                mapping["F'"] = sse["label"]
            if sse["start"] == 201 and sse["end"] == 207: #H20
                mapping["F''"] = sse["label"]
            if sse["start"] == 208 and sse["end"] == 235: #H22
                mapping["G"] = sse["label"]
            if sse["start"] == 244 and sse["end"] == 258:
                mapping["H"] = sse["label"]
            if sse["start"] == 265 and sse["end"] == 297:
                mapping["I"] = sse["label"]
            if sse["start"] == 298 and sse["end"] == 312:
                mapping["J"] = sse["label"]
            if sse["start"] == 327 and sse["end"] == 340:
                mapping["K"] = sse["label"]
            if sse["start"] == 349 and sse["end"] == 350:
                mapping["1_4"] = sse["label"]
            if sse["start"] == 372 and sse["end"] == 377:
                mapping["K'"] = sse["label"]
            if sse["start"] == 419 and sse["end"] == 436:
                mapping["L"] = sse["label"]
            if sse["start"] == 437 and sse["end"] == 441:#
                mapping["3_3"] = sse["label"]
            if sse["start"] == 466 and sse["end"] == 470:
                mapping["3_2"] = sse["label"]
#            sses_list[sse["label"]] = {"start": sse["start_vector"], "end": sse["end_vector"]}
    return mapping


def sses_model_mapping_1akd(model_annotation_file):
    # maping based on  https://sestra.ncbr.muni.cz/job/cyp-1akd
    mapping = {}
    with open(model_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = f"{MODEL[0]}{MODEL[1]}{MODEL[2]}"
        for sse in data[pcd]["secondary_structure_elements"]:
            if sse["start"] == 34 and sse["end"] == 36:
                mapping["A'"] = sse["label"]
            if sse["start"] == 38 and sse["end"] == 46:
                mapping["A"] = sse["label"]
            if sse["start"] == 53 and sse["end"] == 56:
                mapping["1_1"] = sse["label"]
            if sse["start"] == 68 and sse["end"] == 77:
                mapping["B"] = sse["label"]
            if sse["start"] == 90 and sse["end"] == 95:
                mapping["B'"] = sse["label"]
            if sse["start"] == 105 and sse["end"] == 119:
                mapping["C"] = sse["label"]
            if sse["start"] == 121 and sse["end"] == 145:
                mapping["D"] = sse["label"]
            if sse["start"] == 150 and sse["end"] == 167:
                mapping["E"] = sse["label"]
            if sse["start"] == 171 and sse["end"] == 186:
                mapping["F"] = sse["label"]
            if sse["start"] == 193 and sse["end"] == 214:
                mapping["G"] = sse["label"]
            if sse["start"] == 219 and sse["end"] == 225:
                mapping["H"] = sse["label"]
            if sse["start"] == 235 and sse["end"] == 267:
                mapping["I"] = sse["label"]
            if sse["start"] == 268 and sse["end"] == 277:
                mapping["J"] = sse["label"]
            if sse["start"] == 278 and sse["end"] == 292:
                mapping["K"] = sse["label"]
            if sse["start"] == 298 and sse["end"] == 301:
                mapping["1_4"] = sse["label"]
            if sse["start"] == 322 and sse["end"] == 327:
                mapping["K'"] = sse["label"]
            if sse["start"] == 360 and sse["end"] == 378:
                mapping["L"] = sse["label"]
            if sse["start"] == 382 and sse["end"] == 383:
                mapping["3_3"] = sse["label"]
            if sse["start"] == 403 and sse["end"] == 405:
                mapping["3_2"] = sse["label"]
#            sses_list[sse["label"]] = {"start": sse["start_vector"], "end": sse["end_vector"]}
    return mapping


def test_domain_completeness(domain_annotation_file, mapping, percentage=80):
    with open(domain_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = list(data.keys())[0]
        count_found = 0
        found_list = []
        wanted_list = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
        for wanted in  wanted_list:
            for sse in data[pcd]["secondary_structure_elements"]:
                if sse["label"] == mapping[wanted]:
                    count_found += 1
                    found_list.append(wanted)
    print(f"count_found  {count_found} from {len(wanted_list)}, ({found_list})")
   
    if count_found > len(wanted_list)*percentage/100:
        return True
    else:
        return False

def rotate_and_translate(A: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    ''' Applie rotation (or improper rotation) matrix 'rotation' and translation vector 'translation'
    to a 2D object A,
    where columns of A are coordinates of individual points.
    A_result = R * A + t
    '''
    if isinstance(A, dict):
        # on dictionaries:
        matrix = np.array(list(A.values())).transpose()
        matrix = rotate_and_translate(matrix, rotation, translation)
        return { name: list(matrix[:,i]) for i, name in enumerate(A.keys()) }
    else:
        # on matrices:
        if A.size == 0:
            return A
        else:
            return np.matmul(rotation, A) + translation


def predict_hem_position(model_annotation, prot_annotation, cif_filename, rotation_json, model_rotation, model_cif):

    # at first find the annotation
    pos = sses_model_mapping_1og2(model_annotation)

    # find out the position of relevant SSEs which define HEM position
    with open(prot_annotation, "r") as fp:
        prot_data = json.load(fp)
        prot_name = list(prot_data.keys())[0]

        # we need helices B, I and L
        if "I" not in pos or "L" not in pos or "B" not in pos:
            return None

        coord = {}
        res = {}
        for sse in prot_data[prot_name]["secondary_structure_elements"]:
            if (sse["label"] == pos["I"]):
                coord["I"] = [(sse["start_vector"][0]+sse["end_vector"][0])/2, \
                              (sse["start_vector"][1]+sse["end_vector"][1])/2, \
                              (sse["start_vector"][2]+sse["end_vector"][2])/2]
                print("start I", prot_annotation, sse["start_vector"])
                res["I"] = int((sse['start']+sse['end'])/2)
                print(".end. I", prot_annotation, sse["end_vector"])
            if (sse["label"] == pos["L"]):
                coord["L"] = [sse["start_vector"][0], sse["start_vector"][1], sse["start_vector"][2]]
                res["L"] = sse["start"]
            if (sse["label"] == pos["B"]):
                coord["B"] = [sse["end_vector"][0], sse["end_vector"][1], sse["end_vector"][2]]

#    return [(coord["I"][0]+coord["L"][0]+coord["B"][0])/3, (coord["I"][1]+coord["L"][1]+coord["B"][1])/3, (coord["I"][2]+coord["L"][2]+coord["B"][2])/3]
#    print("res",res["I"], utils_read_cif.find_res_avg_coord(cif_filename, "A", res["I"]))
    cc = {}
    print("coord - stred I", prot_annotation, coord["I"], res["I"])
    cc["I"] = utils_read_cif.find_res_avg_coord(cif_filename, "A", res["I"])
    print("coord - start L", prot_annotation, coord["L"], res["L"])
    cc["L"] = utils_read_cif.find_res_avg_coord(cif_filename, "A", res["L"])
    print("coord - stred I", prot_annotation, coord["I"])
    print("coord - start L", prot_annotation, coord["L"])
    print("old vysledek", prot_annotation, [(coord["I"][0]+coord["L"][0])/2, (coord["I"][1]+coord["L"][1])/2, (coord["I"][2]+coord["L"][2])/2])
    model_hem = utils_read_cif.find_fe_hem_coords_auth(model_cif, MODEL[1])

    with open(model_rotation, "r") as fr:
        r_data = json.load(fr)
#        m_r = np.linalg.inv(r_data["rotation"])
        m_r = np.array(r_data["rotation"])
        m_t = np.array(r_data["translation"])
#        m_r2 = np.linalg.inv(r_data["rotation"])

    with open(rotation_json, "r") as fr:
        r_data = json.load(fr)
        r_inv = np.linalg.inv(r_data["rotation"])
        r = np.array(r_data["rotation"])
        t = np.array(r_data["translation"])
#        tt = np.array([t[0][0], t[1][0], t[2][0]])
#        hem = np.array([(coord["I"][0]+coord["L"][0])/2, (coord["I"][1]+coord["L"][1])/2, (coord["I"][2]+coord["L"][2])/2])
#        hem_t = [hem[0]-t[0][0], hem[1] - t[1][0], hem[2]-t[2][0]]
#        transf_hem = [r[0][0]*hem[0]+r[1][0]*hem[1]+r[2][0]*hem[2]+t[0][0], \
#                     r[0][1]*hem[0]+r[1][1]*hem[1]+r[2][1]*hem[2]+t[1][0], \
#                     r[0][2]*hem[0]+r[1][2]*hem[1]+r[2][2]*hem[2]+t[2][0]]
    print("------")
    print("vysledek", cif_filename)
#    print("r", type(r), r)
#    print("hem", type(hem), hem)
#    aux = np.dot(r, hem)
#    print("dot", type(aux), aux)
    
#    print("tt", type(tt), tt)
#    print("22222222",  np.add(aux, tt))
#    print("start")
    print("model hem po rotaci 2:", np.dot(m_r.T, np.array(model_hem).T))
#    print("m_t translace", m_t, m_t.T)
    print("po rotaci translaci 2", np.dot(m_r, np.array(model_hem)) + m_t.T)
#    print("po rotaci translaci 3", (np.dot(m_r2.T, np.array(model_hem).T) + m_t.T))
    transf_model_hem = np.dot(m_r, np.array(model_hem)) + m_t.T
    print("transf_model_hem", transf_model_hem)
    print("t", t)
    print("r_inv", r_inv)
#    print("11111111", transf_model_hem.T-t)
    ret = np.dot(np.linalg.inv(r), transf_model_hem.T-t)
    rfl = list(ret.flatten())
    print("33333333", ret, [float(rfl[0]), float(rfl[1]), float(rfl[2])])
#    print("44444444", ret, list((ret.T)[0]),)
#    print("55545555", ret, [ret[0], ret[1], ret[2]])
    return [float(rfl[0]), float(rfl[1]), float(rfl[2])]
#    sys.exit()
#    sys.exit()
        
    
#    return [(coord["I"][0]+coord["L"][0])/2, (coord["I"][1]+coord["L"][1])/2, (coord["I"][2]+coord["L"][2])/2]
    return cc["I"]
#    return [(cc["I"][0]+cc["L"][0])/2, (cc["I"][1]+cc["L"][1])/2, (cc["I"][2]+cc["L"][2])/2]
