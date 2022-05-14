import requests
import xml.etree.ElementTree as ET
import time
import concurrent.futures
from tqdm import tqdm
import pandas as pd

ENTREZ_SEARCH_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

API_KEY = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"  # Note: this key has been disabled

params_base = {
    'db': 'pubmed',
    'tool': 'GeneID',
    'email': 'david.detomaso@berkeley.edu',
    'api_key': API_KEY,
}

scope = "[Title/Abstract]"


def mk_term_gene(gene, context=None):
    """
    context is a list of terms to include with gene
    search is gene AND (context[0] OR context[1] OR ...)
    """

    g_string = gene + scope

    if context is not None:
        c_string = " OR ".join(['"{}"{}'.format(x, scope) for x in context])
        c_string = "({})".format(c_string)
        res = g_string + " AND " + c_string
    else:
        res = g_string

    return res


def gene_count(gene, context=None):

    query = mk_term_gene(gene, context)

    params = params_base.copy()
    params['term'] = query

    r = requests.get(ENTREZ_SEARCH_BASE, params=params)

    # Get the count from the result
    root = ET.fromstring(r.text)
    c = int(root.find('Count').text)
    return c


def gene_count_rl(gene, rate_limit, context=None):
    t1 = time.time()

    try:
        c = gene_count(gene, context)
    except Exception:
        c = -1

    t2 = time.time()
    if t2 - t1 < rate_limit:
        time.sleep(rate_limit - (t2 - t1))

    return (gene, c)


def get_gene_counts(genes, context=None):
    N_WORKERS = 10
    RATE_LIMIT = 1.5  # limit in seconds for the N worksers

    res = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=N_WORKERS) as executor:
        results = [
            executor.submit(gene_count_rl, gene, RATE_LIMIT, context) for gene in genes]
        for future in tqdm(concurrent.futures.as_completed(results), total=len(genes)):
            gene, count = future.result()
            res[gene] = count

    res = pd.Series(res, name='Count')

    return res


tpm_in = snakemake.input['tpm']
out = snakemake.output['out']

# Load data
tpm = pd.read_table(tpm_in, index_col=0)
genes = tpm.index.tolist()

c_base = get_gene_counts(genes)

context = [
    "Th1",
    "CD4",
    "T Cell",
    "lymphocyte",
]

c_tcell = get_gene_counts(genes, context)

res = pd.concat(
    (c_base.rename('Base'), c_tcell.rename('TCell')),
    axis=1, sort=False
)

res.to_csv(out, sep="\t")
