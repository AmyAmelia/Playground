from pybiomart import Server, Dataset

#
server = Server(host='http://www.grch37.ensembl.org')
# server = Server(host='http://www.grch37.ensembl.org')
print(server)

mart = server['ENSEMBL_MART_ENSEMBL']
print(mart.list_datasets())


# dataset = Dataset(name='hsapiens_gene_ensembl')
# print(dataset.query())
#
# print(dataset)
#
# dataset = (server.marts['ENSEMBL_MART_ENSEMBL']
#             .datasets['hsapiens_gene_ensembl'])

#
# dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
#               filters={'chromosome_name': ['1','2']})


# dataset = Dataset(name='hsapiens_gene_ensembl', host = 'http://www.ensembl.org')
# dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
#
# dataset.attributes
