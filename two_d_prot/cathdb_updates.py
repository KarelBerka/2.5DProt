import logging

import sqlalchemy
import sqlalchemy.orm
import sqlalchemy.ext.declarative

import two_d_prot.cathdb


Base = sqlalchemy.ext.declarative.declarative_base()


class Family(Base):
    __tablename__ = 'families'

    name = sqlalchemy.Column(sqlalchemy.String(50), primary_key=True)
    last_updated_at = sqlalchemy.Column(sqlalchemy.Date())
    last_deleted_at = sqlalchemy.Column(sqlalchemy.Date())

    def __repr__(self):
        if self.last_deleted_at is None:
            return "<Family(name='%s', last_updated_at='%s')>" % (
                self.name, self.last_updated_at)
        else:
            return "<Family(name='%s', last_deleted_at='%s')>" % (
                self.name, self.last_deleted_at)


class Domain(Base):
    __tablename__ = 'domains'

    name = sqlalchemy.Column(sqlalchemy.String(20), primary_key=True)
    family = sqlalchemy.Column(sqlalchemy.String(50), sqlalchemy.ForeignKey('families.name'))
    residues = sqlalchemy.Column(sqlalchemy.String(200))
    last_updated_at = sqlalchemy.Column(sqlalchemy.Date())

    def __repr__(self):
        return "<Domain(name='%s', family='%s', residues='%s', last_updated_at='%s')>" % (
            self.name, self.family, self.residues, self.last_updated_at)


class CathUpdates(object):

    def __init__(self):
        SessionMaker = self.open_or_create_db()
        self.session = SessionMaker()

    def open_or_create_db(self):
        engine = sqlalchemy.create_engine('sqlite:///input/cath-family-updates.sqlite', echo=False)
        engine.connect()

        Base.metadata.create_all(engine)

        SessionMaker = sqlalchemy.orm.sessionmaker(bind=engine)

        return SessionMaker

    def show_info(self):
        for family in self.session.query(Family).order_by(Family.name):
            print("%s" % family)
            for domain in self.session.query(Domain).filter(Domain.family == family.name).order_by(Domain.name):
                print("    %s" % domain)

    def last_updated_per_family(self):
        families = self.session.query(Family) \
            .filter(Family.last_deleted_at == None) \
            .order_by(Family.last_updated_at).order_by(Family.name)   # noqa: E711
        for family in families:
            yield (family.name, family.last_updated_at)

    def family_updated(self, family):
        family_instance = self.session.query(Family).filter(Family.name == family).first()
        return family_instance.last_updated_at.isoformat()

    def families_after(self, after):
        families = self.session.query(Family) \
            .filter(sqlalchemy.and_(Family.last_deleted_at == None, Family.last_updated_at >= after)) \
            .order_by(Family.name)   # noqa: E711
        for family in families:
            yield family.name

    def domains_after(self, after, family):
        if family is None:
            domains = self.session.query(Domain).filter(Domain.last_updated_at >= after).order_by(Domain.name)
        else:
            domains = self.session.query(Domain).filter(Domain.family == family).filter(Domain.last_updated_at >= after).order_by(Domain.name)
        for domain in domains:
            yield (domain.family, domain.name)

    def families_size(self):
        families = self.session.query(Family) \
            .filter(Family.last_deleted_at == None) \
            .order_by(Family.name)   # noqa: E711
        for family in families:
            count = self.session.query(Domain).filter(Domain.family == family.name).count()
            yield (family.name, count, family.last_updated_at)

    def last_db_update(self):
        return self.session.query(sqlalchemy.func.max(Family.last_updated_at)).first()[0].strftime("%Y-%m-%d")

    def create_domain(self, domain, family, residues, last_updated_at):
        domain_instance = self.session.query(Domain).filter(Domain.name == domain).first()
        if domain_instance is None:
            domain_instance = Domain(name=domain, family=family, residues=residues,
                                     last_updated_at=last_updated_at)
            self.session.add(domain_instance)
            logging.debug("For family %s creating %s" % (family, domain_instance))
        else:
            logging.debug("Moving domain from family %s to family %s: %s" % (domain_instance.family, family, domain_instance))
            family_instance = self.session.query(Family).filter(Family.name == domain_instance.family)
            family_instance.last_updated_at = last_updated_at
            domain_instance.family = family

    def update_from_release(self, date, putative_only):
        url = "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/daily-release/archive/cath-b-%s-all.gz" % date.strftime('%Y%m%d')
        cath = two_d_prot.cathdb.CathDailyRelease(url)
        if putative_only:
            families = cath.families_putative
        else:
            families = cath.families
        logging.info("Working with %s families in putative mode? %s" % (len(families), putative_only))
        counter = 0
        for family in families:
            counter += 1
            if counter % 100 == 0:
                pending = len(self.session.new) + len(self.session.dirty) + len(self.session.deleted)
                logging.debug("Detected %s pending objects in session" % pending)
                if pending > 1000:
                    self.session.commit()
                    logging.debug("Session commited")

            current_instance = self.session.query(Family).filter(Family.name == family).first()

            new_residues = {d.name: d.boundaries for d in cath.get_domains_for_family(family)}
            new_updated_at = date

            # If we are just adding new family into DB:
            if current_instance is None:
                family_instance = Family(
                    name=family, last_updated_at=new_updated_at,
                    last_deleted_at=None)
                self.session.add(family_instance)
                logging.debug("Creating new %s" % family_instance)
                for k, v in new_residues.items():
                    self.create_domain(k, family, v, new_updated_at)
                continue

            current_updated_at = current_instance.last_updated_at

            # If update we are trying to push is older than what we already have:
            if new_updated_at < current_updated_at:
                logging.debug("Too old change for family %s" % family)
                continue

            # Load current residues
            current_residues = {}
            for d in self.session.query(Domain).filter(Domain.family == family).order_by(Domain.name).all():
                current_residues[d.name] = d.residues
            logging.debug("Loaded %s domains for family %s" % (len(current_residues), family))

            # Skipp it if after all there was no change, othervise we know we are
            # going to update the family so update its last_updated_at
            if current_residues == new_residues:
                logging.debug("No change in domains and residues for family %s" % family)
                continue
            else:
                logging.debug("Setting new last updated date %s for family %s" % (new_updated_at, family))
                current_instance.last_updated_at = new_updated_at

                # So, if something really changed

                # ... new domains
                for d in set(new_residues.keys()) - set(current_residues.keys()):
                    self.create_domain(d, family, new_residues[d], new_updated_at)

                # ... changed domains
                for d in set(current_residues.keys()).intersection(set(new_residues.keys())):
                    if current_residues[d] != new_residues[d]:
                        domain_instance = self.session.query(Domain).filter(Domain.name == d).first()
                        domain_instance.residues = new_residues[d]
                        domain_instance.last_updated_at = new_updated_at
                        logging.debug("For family %s updated residues for %s" % (family, domain_instance))

                # ... removed domains
                for d in set(current_residues.keys()) - set(new_residues.keys()):
                    domain_instance = self.session.query(Domain).filter(Domain.name == d).first()
                    self.session.delete(domain_instance)
                    logging.debug("For family %s deleting %s" % (family, domain_instance))

        self.session.commit()

        logging.info("Checking all (not only putative) families if these need to be deleted")

        new_families = cath.families
        current_families = self.session.query(Family) \
            .filter(Family.last_deleted_at == None)
        current_families_names = [f.name for f in current_families]

        for family in set(current_families_names) - set(new_families):
            # Find family instance
            for family_instance in current_families:
                if family_instance.name == family:
                    break

            # Remove its domains and family itself
            for domain_instance in self.session.query(Domain).filter(Domain.family == family).all():
                self.session.delete(domain_instance)
            self.session.delete(family_instance)
            logging.debug("Deleted family %s" % family_instance)

            self.session.commit()
