import logging
import json
import azure.functions as func


def main(req: func.HttpRequest, documents: func.DocumentList) -> func.HttpResponse:
    logging.info('Python HTTP trigger function processed a request.')

    humanCount = None
    mutantCount = None
    
    for countDocument in documents:
        countDict = json.loads(countDocument.to_json())
        if countDict['DnaType'] == 'mutant':
            mutantCount = countDict['Count']
        elif countDict['DnaType'] == 'human':
            humanCount = countDict['Count']

    logging.info(f"from database:{{\"count_mutant_dna:{mutantCount}\",\"count_human_dna\":{humanCount}")

    # if count from query are OK
    if humanCount == None: humanCount = 0
    if mutantCount == None: mutantCount = 0

    if humanCount == 0:
        if mutantCount > 0:
            ratio = mutantCount/mutantCount
        else: ratio = 0
    else:
        ratio = mutantCount/humanCount

    response = f"{{\"count_mutant_dna\":{mutantCount},\"count_human_dna\":{humanCount},\"ratio\":{ratio}}}"

    logging.info(f"Response string: {response}")
    
    return func.HttpResponse(response)
