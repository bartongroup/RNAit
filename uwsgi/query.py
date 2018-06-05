#!/usr/bin/env python

from cgi import parse_qs
from Bio import SeqIO
import io
import primer3
import pprint

import cgitb
cgitb.enable()

def application (environ,start_response):
    start_response('200 OK',[('Content-Type','text/html')])
    
    params = get_params(environ)
    if ('error' in params):
        html = params['error']
        
        encode=html.encode('UTF-8')
        return(encode)
    
    primers = get_primer_pairs(params)
    
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