import logging
import numpy as np
import itertools as it
import json
import azure.functions as func

def main(req: func.HttpRequest, outputDocument: func.Out[func.Document]) -> func.HttpResponse:
    logging.info('Python HTTP trigger function processed a request.')
    
    dna = req.params.get('dna')
    if not dna:
        try:
            req_body = req.get_json()
        except ValueError:
            pass
        else:
            dna = req_body.get('dna')

    if dna:
        try:
            dnaData = json.loads(dna)
            dnaData = [value for value in dnaData]
            dnaData = [list(s) for s in dnaData]
        except:
            return func.HttpResponse(
             "Please pass mutant dna on the query string or in the request body",
             status_code=400
            )

        try:
            matches = mutant_count(np.array(dnaData))
        except:
            return func.HttpResponse(
             "Please pass mutant dna on the query string or in the request body",
             status_code=500
            )
       
        logging.info(f'DNA:{dnaData} match count:{matches}.')

        if (matches == None):
            return func.HttpResponse(
                "Please pass mutant dna on the query string or in the request body",
                status_code=400
            )
        
        if (matches > 1):
            pk = 'mutant'
            status_code=200
        elif (matches == 0):
            pk = 'human'
            status_code=403

        data = {
                "PartitionKey": pk,
                "id": dna,
                "RowKey": dna,
                "dna": dna
            }

        outputDocument.set(func.Document.from_dict(data))
        return func.HttpResponse(
                status_code=status_code
            )
    else:
        return func.HttpResponse(
             "Please pass mutant dna on the query string or in the request body",
             status_code=400
        )


def mutant_count(adn_matrix, adn_base=['A', 'C', 'G', 'T'], adn_seq_length=4):
    N = len(adn_matrix)

    # verify if DNA sequence contains valid DNA elements
    ACGT_count = np.count_nonzero(np.logical_or(adn_matrix == 'A', adn_matrix == 'T')) +  np.count_nonzero(np.logical_or(adn_matrix ==  'C', adn_matrix ==  'G'))

    if (ACGT_count != N*N):
        return None

    match_count = 0

    # l: letter, g: group
    # main diagonals
    match_count += np.sum(
        [len(list(g))  >= adn_seq_length 
        for l, g in it.groupby(np.diagonal(adn_matrix,offset=0)) 
        if l in adn_base])

    match_count += np.sum(
        [len(list(g))  >= adn_seq_length 
        for l, g in it.groupby(np.fliplr(adn_matrix).diagonal(offset=0)) 
        if l in adn_base])

    # Diagonals without main diagonals
    max_diag_size = range(1, N-adn_seq_length+1 )
    for i in range(1, N-adn_seq_length+1): 
        # l: letter, g: group
        # left to right upper diagonals
        match_count += np.sum(
            [len(list(g))  >= adn_seq_length 
            for l, g in it.groupby(np.diagonal(adn_matrix,offset=i)) 
            if l in adn_base])
        
        # # left to right lower diagonals
        match_count += np.sum(
            [len(list(g))  >= adn_seq_length 
            for l, g in it.groupby(np.diagonal(adn_matrix,offset=-i)) 
            if l in adn_base])

        #right to left upper diagonals
        match_count += np.sum(
            [len(list(g))  >= adn_seq_length 
            for l, g in it.groupby(np.fliplr(adn_matrix).diagonal(offset=i)) 
            if l in adn_base])

        # #right to left lower diagonals
        match_count += np.sum(
            [len(list(g))  >= adn_seq_length 
            for l, g in it.groupby(np.fliplr(adn_matrix).diagonal(offset=-i)) 
            if l in adn_base])
        
    # Remaining rows and columns
    for i in range(1, N+1):
        match_count += np.sum(
            [len(list(g))  >= adn_seq_length 
            for l, g in it.groupby(adn_matrix[i-1,:]) 
            if l in adn_base])

        match_count += np.sum(
            [len(list(g))  >= adn_seq_length 
            for l, g in it.groupby(adn_matrix[:,i-1]) 
            if l in adn_base])
        
    return int(match_count)