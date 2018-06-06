#!/usr/bin/env python

from cgi import parse_qs
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from tempfile import mkdtemp
import io
import re
import primer3
import pprint
import shutil

import cgitb
cgitb.enable()

def application (environ,start_response):
    start_response('200 OK',[('Content-Type','text/html')])
    pp=pprint.PrettyPrinter(indent=4)
    
    params = get_params(environ)
    if ('error' in params):
        html = params['error']
        
        encode=html.encode('UTF-8')
        return(encode)
    
    seq = params.get('seq')
    db = params.get('database')
    primers = get_primer_pairs(params)
    
    if (len(primers)==0):
        html = 'No suitable primers found'
        encode=html.encode('UTF-8')
        return(encode)
    for pair in primers:
        product=get_pcr_product(seq,pair)
        blast_results=blast_product(product,db)
        
    html=''
    encode=html.encode('UTF-8')
        
    return [encode]


# get_params
#
# parses form parameters following form submission
# Sequences entered via the 'seqpaste' field are parsed
# as SeqRecord objects via SeqIO
#
# required args: environ - environment dictionary
#
# returns: params - dictionary of parsed parameters
            
# TODO: validation of submitted values

def get_params(environ):
    
    try:
        request_body_size = int(environ.get('CONTENT_LENGTH',0))
    except (ValueError):
        request_body_size = 0
    request_body = environ['wsgi.input'].read(request_body_size) 
    d = parse_qs(request_body)
    
    params={}
    for key in d.keys():
        val = d.get(key,[''])[0].decode("utf-8")
        if (key.decode("utf-8")=='seqpaste'):
            seqH = io.StringIO(val)
            try:
                record = SeqIO.read(seqH,'fasta')
                params['seq']=record
            except ValueError:
                params['error']='The entered sequence does not appear to be valid fasta format'
        else:
            params[key.decode("utf-8")]=val
        
    return(params)

# get_primer_pairs
#
# runs primer3 using primer3-py bindings against the provided sequence according to the requested
# melting temperature and product size parameters
#
# required args: params - dictionary of parsed form parameters
#
# returns: primers - list of primer pair dictionaries

def get_primer_pairs(params):
    
    seq_args={
        'SEQUENCE_ID': params.get('seq').id,
        'SEQUENCE_TEMPLATE': str(params.get('seq').seq),
    }
    
    global_args={
        'PRIMER_TASK': 'generic',
        'PRIMER_PRODUCT_SIZE_RANGE': [[int(params.get('product_min')),int(params.get('product_max'))]],
        'PRIMER_OPT_TM': int(params.get('melting_temp')),
    }
    
    primers=primer3.bindings.designPrimers(seq_args,global_args)
    
    pair_count = int(primers.get('PRIMER_PAIR_NUM_RETURNED'))
    pairs=[]
    for i in range(pair_count):
        pair={
            'LEFT_START':(primers.get('PRIMER_LEFT_'+str(i)))[0],
            'LEFT_LENGTH':(primers.get('PRIMER_LEFT_'+str(i)))[1],
            'RIGHT_START':(primers.get('PRIMER_RIGHT_'+str(i)))[0],
            'RIGHT_LENGTH':(primers.get('PRIMER_RIGHT_'+str(i)))[1],
            'LEFT_SEQ': primers.get('PRIMER_LEFT_'+str(i)+'_SEQUENCE'),
            'RIGHT_SEQ': primers.get('PRIMER_RIGHT_'+str(i)+'_SEQUENCE'),
            'LEFT_GC': primers.get('PRIMER_LEFT_'+str(i)+'_GC_PERCENT'),
            'RIGHT_GC': primers.get('PRIMER_RIGHT_'+str(i)+'_GC_PERCENT'),
            'LEFT_MELTING': primers.get('PRIMER_LEFT_'+str(i)+'_TM'),
            'RIGHT_MELTING': primers.get('PRIMER_RIGHT_'+str(i)+'_TM'),
            'LEFT_END_STAB': primers.get('PRIMER_LEFT_'+str(i)+'_END_STABILITY'),
            'RIGHT_END_STAB': primers.get('PRIMER_RIGHT_'+str(i)+'_END_STABILITY'),
            'PRODUCT_SIZE': primers.get('PRIMER_PAIR_'+str(i)+'_PRODUCT_SIZE'),
            'COMP_END': primers.get('PRIMER_PAIR_'+str(i)+'_COMPL_END'),
        }
        pairs.append(pair)
        
    return(pairs)

# get_pcr_product
#
# Isolates the subsequence represnting the pcr product for a primer pair
#
# required args: seq - Bio.seqRecord object
#                pair - dictionary of primer pair info
#
# returns: product -  Bio.seqRecord object

def get_pcr_product(seq,pair):
    
    #both primer3 and biopython use 0 based co-ordinates so we should be good to go....
    start=pair.get('LEFT_START')
    end=pair.get('RIGHT_START')+1
    product=seq[start:end]
    return(product)

# blast_product
#
# Blasts pcr product against organism genome database to identify
#
# required args: product - Bio:seqRecord object representing pcr product
#                db - blast database name
#
# returns: blast_data - dictionary containing 'record' (blast_record object),
#           alignment_status (status of each alignment [same,conflicting,ok])
#           status ([bad/suitable])

def blast_product(product,db):
    tmpdir=mkdtemp(dir='../tmp')
    queryFileName=tmpdir+'/query'
    outFileName=tmpdir+'/output.xml'
    
    SeqIO.write(product,queryFileName,'fasta')
    cline = NcbiblastnCommandline(cmd='blastn', query=queryFileName, out=outFileName, outfmt=5, db='../databases/'+db, evalue=0.01)
    stdout,stderr=cline()
    
    result_handle = open(outFileName)
    status=''
    
    blast_record = NCBIXML.read(result_handle)
    midline_regex=re.compile(r"\|{20,}") 
    alignment_status=[]
    for alignment in blast_record.alignments:
        min_ident=1.00
        status_set=0
        conflicting=0
        for hsp in alignment.hsps:
            # check for matches of >20bp identity by checking for stretches of >20 '|' characters in the HSP midline
            match=midline_regex.search(hsp.match)
            match_len=match.end()-match.start()
            ident=hsp.identities/hsp.align_length
            if (ident>0.99):
                alignment_status.append('same')
                status_set=1
            if (min_ident>ident and (ident>0.89 and ident<0.99)):
                min_ident=ident
            
        if ((min_ident>0.89 or match_len>20) and status_set==0):
            alignment_status.append('conflicting: id='+str(min_ident*100)+'; identical stretch='+str(match_len))
            conflicting=1
        elif (status_set==0):
            aligment_status.append('ok')
            
    if conflicting:
        status='bad'
    else:
        status='suitable'
    
    blast_data={
        'record': blast_record,
        'alignment_status': alignment_status,
        'primer_status': status
        }
    shutil.rmtree(tmpdir)
    
    return(blast_data)
    



