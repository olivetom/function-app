import logging
import numpy as np
import itertools as it
import json
import azure.functions as func

def main(req: func.HttpRequest) -> func.HttpResponse:
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
        response = json.loads(dna)
        response = [value for value in response]
        response = [list(s) for s in response]
        matches = mutant_count(np.array(response))
        body=f"match_count:{matches}"
        if (matches > 1):
            return func.HttpResponse(
                body=body,
                status_code=200
            )
        else:
            return func.HttpResponse(
                body=body,
                status_code=403
            )
    else:
        return func.HttpResponse(
             "Please pass mutant dna on the query string or in the request body",
             status_code=400
        )


def mutant_count(adn_matrix, adn_base=['A', 'C', 'G', 'T'], adn_seq_length=4):
    N = len(adn_matrix)
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