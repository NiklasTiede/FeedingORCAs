
import pymongo

# MongoDB molecule collection used within feedingORCAs:
from typing import List

DB_NAME = "test_database"
COLLECTION_NAME = "molecules"
# DB_NAME = "feedingORCAs_MongoDB"
# COLLECTION_NAME = "Molecules"


class MoleculeMongoDB:

    def __init__(self):
        # connects to default host/port:
        self.client = pymongo.MongoClient()
        # select a db and collection (eq. to sql-table)
        self.db = self.client[DB_NAME]
        self.molecules_doc = self.db[COLLECTION_NAME]

    def query_all(self) -> List[dict]:
        documents = []
        for doc in self.molecules_doc:
            documents.append(doc)
        return documents


# from ..molecular.molecule import Molecule




