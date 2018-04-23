"""
Get set of edge types in the graph
"""
from neo4j.v1 import GraphDatabase

query = "call db.relationshipTypes"
rel_types = []

scigraph_bolt = "bolt://neo4j.monarchinitiative.org:443"
driver = GraphDatabase.driver(scigraph_bolt, auth=("neo4j", "neo4j"))

def get_scigraph_results(query):
    with driver.session() as session:
        with session.begin_transaction() as tx:
            for record in tx.run(query):
                yield record

for record in get_scigraph_results(query):
    rel_types.append(record['relationshipType'])

for rel in rel_types:
    label_query = """
    MATCH ()-[r:`{}`]-()
    RETURN r.lbl as label
    LIMIT 1
    """.format(rel)
    for lbl in get_scigraph_results(label_query):
        print("{}\t{}".format(rel, lbl['label']))