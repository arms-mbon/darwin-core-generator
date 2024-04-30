from argparse import Namespace
import ncbitaxonomist.groupmanager
import ncbitaxonomist.log.logger
import ncbitaxonomist.mapper
import ncbitaxonomist.collector
import ncbitaxonomist.ncbitaxonomist
import ncbitaxonomist.parser.arguments
import ncbitaxonomist.payload.accession
import ncbitaxonomist.payload.name
import ncbitaxonomist.payload.taxid
import ncbitaxonomist.resolve.resolver
import ncbitaxonomist.subtree.subtreeanalyzer
import ncbitaxonomist.utils

def resolve(name):
    args = Namespace()
    args.version = False
    args.verbose = 0
    args.apikey = None
    args.command = "resolve"
    args.taxids = None
    args.names = [name]
    args.database = None
    args.remote = True
    args.email = None
    args.xml = False
    args.mapping = False

    try:
        ncbitaxonomist.ncbitaxonomist.configure(args)
        nt = ncbitaxonomist.ncbitaxonomist.NcbiTaxonomist(args.database)
        txresolver = ncbitaxonomist.resolve.resolver.Resolver(nt)
        txresolver.cache.taxa.taxa = {}
        txresolver.resolve(
            taxids=ncbitaxonomist.payload.taxid.TaxidPayload(args.taxids),
            names=ncbitaxonomist.payload.name.NamePayload(args.names),
            mapping=args.mapping,
            remote=args.remote,
        )
        ncbi_id = [k for k, _ in txresolver.cache.taxa.taxa.items()][0]
        return ncbi_id
    except Exception:
        return None
