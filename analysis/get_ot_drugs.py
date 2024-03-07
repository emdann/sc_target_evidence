import requests
import pandas as pd

def download_drug_target_data(efo_id):
    url = 'https://api.platform.opentargets.org/api/v4/graphql'

    headers = {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
    }

    query = """
    query KnownDrugsQuery(
      $efoId: String!
      $cursor: String
      $freeTextQuery: String
      $size: Int = 10
    ) {
      disease(efoId: $efoId) {
        id
        knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
          count
          cursor
          rows {
            phase
            status
            urls {
              name
              url
            }
            disease {
              id
              name
            }
            drug {
              id
              name
              mechanismsOfAction {
                rows {
                  actionType
                  targets {
                    id
                  }
                }
              }
            }
            urls {
              url
              name
            }
            drugType
            mechanismOfAction
            target {
              id
              approvedName
              approvedSymbol
            }
          }
        }
      }
    }
    """

    variables = {
        "efoId": efo_id
    }

    response = requests.post(url, json={'query': query, 'variables': variables}, headers=headers)

    if response.status_code == 200:
        data = response.json()['data']['disease']['knownDrugs']['rows']
        df = pd.json_normalize(data)
        return df
    else:
        print(f"Error: Unable to fetch data. Status code: {response.status_code}")
        return None

# Example usage:
efo_id = "EFO:0000305"  # Example EFO ID for a disease (e.g., asthma)
df = download_drug_target_data(efo_id)
print(df)
