#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Example:
from search_rcsb import searchMacromolecule

searchMacromolecule('adrenergic receptor')
searchMacromolecule('rhodopsin')
'''

import requests
from urllib.parse import quote

def parseRequest(query, viewOnline = False):
    # The base url to request a query...
    BASE = "https://search.rcsb.org/rcsbsearch/v1/query"
    BASE_replace = "https://www.rcsb.org/search?request="

    # Form the url...
    # Convert to double quote (%22) from single quote (%27) 
    # Convert to "False" from "false"
    query_string = quote(f'{query}').replace("%27", "%22").replace("%22false%22", "false")
    url = f'{BASE}?json={query_string}'

    # Request...
    r = requests.get( url = url )

    if viewOnline: 
        # Request...
        print("Copy and paste the following link to your web browser:")
        print("")
        r = requests.get( url = url )

        # Display the URL for visualizing results in web browser...
        print(r.url.replace(f"{BASE}?json=", BASE_replace))
        print("")
        print("done")

    return r




# [[[ Utilities ]]]
def findSimilar(pdb_entry, viewOnline = True):
    # The query to search similar assemblies...
    query = {
      "query": {
            "type": "terminal",
            "service": "structure",
            "parameters": {
              "value": {
                "entry_id": f"{pdb_entry}",
                "assembly_id": "1"
              },
              "operator": "strict_shape_match"
          }
      },
      "request_options": {
            "pager": {
              "start": 0,
              "rows": 20000,
            }
      },
      "return_type": "assembly"
    }

    r = parseRequest(query, viewOnline = viewOnline)
    return r




def findSimilarSequence(sequence, cutoff, viewOnline = False):
    query = {
      "query": {
        "type": "terminal",
        "service": "sequence",
        "parameters": {
          "evalue_cutoff": 1000000,
          "identity_cutoff": cutoff/100.0,
          "target": "pdb_protein_sequence",
          "value": sequence
        }
      },
      "return_type": "entry",
      "request_options": {
        "pager": {
          "start": 0,
          "rows": 10000
        },
        "scoring_strategy": "sequence",
        "sort": [
          {
            "sort_by": "score",
            "direction": "desc"
          }
        ]
      }
    }

    r = parseRequest(query, viewOnline = viewOnline)
    return r




def macromolecule(macromolecule_type, viewOnline = False):
    query = {
      "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "contains_phrase",
          "negation": "false",
          "value": f"{macromolecule_type}",
          "attribute": "rcsb_polymer_entity.pdbx_description"
        }
      },
      "return_type": "entry",
      "request_options": {
        "pager": {
          "start": 0,
          "rows": 1000000
        },
        "scoring_strategy": "combined",
        "sort": [
          {
            "sort_by": "score",
            "direction": "desc"
          }
        ]
      }
    }

    r = parseRequest(query, viewOnline = viewOnline)
    return r
