import requests
import sys
from os import mkdir
from datetime import datetime, timedelta
from lxml import etree

requests.packages.urllib3.disable_warnings()  # I'm using verify=False, so this suppresses the SSL warnings. Not ideal practice....

# bc763bc6-2013-4a7e-8c20-962244299092 => COVID19 Genomics
# 7b85b189-5b1a-4564-b2e2-e51473265d78 => Proportions Archive

# Config  ########
view_names = ['Smallfigure', 'Mapforppt']
regions = ['USA', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
workbook_name = "Variant_Proportions_Plus_Nowcasting_240103" #Default for testing

###################
def get_dates():
    # Fixed reference date
    reference_date = datetime(2039, 12, 31).date()

    # Get today's date
    today = datetime.now().date()

    # Calculate the difference in days from the reference date
    days_difference = (reference_date - today).days

    # Calculate the result using modulo to achieve a biweekly pattern
    result_date = today + timedelta(days=(days_difference % 14)) #default Nowcast
    second_date = result_date - timedelta(weeks=4) #default weighted

    if(today < result_date - timedelta(days=6)):
        print("adjusint dates")
        result_date = result_date - timedelta(weeks=2)
        second_date = second_date - timedelta(weeks=2)
    print(f"'Nowcast' : {result_date}, 'Weighted': {second_date} ")
    
    return { 'Nowcast' : result_date, 'Weighted': second_date }

    

def get_token():
    url = "https://tableau.edav.cdc.gov/api/3.11/auth/signin"
    headers = {'Content-Type': 'application/xml'}
    
    # Load the content of the signin.xml file
    with open("signin.xml", "rb") as file:
        data = file.read()

    # HTTP POST request
    response = requests.post(url, headers=headers, data=data, verify=False)

    # Parse the XML response
    root = etree.fromstring(response.content)
    
    # Find the 'credentials' element using the namespace
    namespace = {'ns': 'http://tableau.com/api'}
    credentials_elem = root.find('.//ns:credentials', namespaces=namespace)

    # Check if the 'credentials' element was found
    if credentials_elem is not None:
        # Extract the token value
        token = credentials_elem.get('token')
        return token
    else:
        # Handle the case where 'credentials' element is not found
        print("Error: 'credentials' element not found in the XML response.", file=sys.stderr)
        return None

# Function to perform HTTP GET with the obtained token and search for a workbook by name
def get_workbook_id(token, workbook_name):
    url = f"https://tableau.edav.cdc.gov/api/3.11/sites/908e3c46-dd30-401c-afe7-f59d3d16e2d0/workbooks?filter=name:eq:{workbook_name}"
    headers = {'X-Tableau-Auth': token}

    # HTTP request
    response = requests.get(url, headers=headers, verify=False)
    # Parse the XML content
    root = etree.fromstring(response.content)
    namespace = {'ns': 'http://tableau.com/api'}
    # Use XPath to find the 'id' attribute for the workbook with the specified name
    result = root.xpath(f'//ns:workbook[@name="{workbook_name}"]/@id', namespaces=namespace)
    
    # Return the luid
    return result[0]


def get_view_id(token, workbook_id, view_name):
    url = f"https://tableau.edav.cdc.gov/api/3.11/sites/908e3c46-dd30-401c-afe7-f59d3d16e2d0/workbooks/{workbook_id}/views"
    headers = {'X-Tableau-Auth': token}

    # HTTP request
    response = requests.get(url, headers=headers, verify=False)
    # Parse the XML content
    root = etree.fromstring(response.content)
    namespace = {'ns': 'http://tableau.com/api'}
    # Use XPath to find the 'id' attribute for the workbook with the specified name
    result = root.xpath(f'//ns:view[@viewUrlName="{view_name}"]/@id', namespaces=namespace)
    # Return the luid
    return result[0] 

def get_image(token, view_id, otherparam, filename):
    url = f"https://tableau.edav.cdc.gov/api/3.11/sites/908e3c46-dd30-401c-afe7-f59d3d16e2d0/views/{view_id}/image?resolution=high&{otherparam}"
    headers = {'X-Tableau-Auth': token} 
    print("getting: ", url)
    response = requests.get(url, headers=headers, verify=False)
    output_path = filename #adjust
    if response.status_code == 200:
        with open(output_path, 'wb') as output_file:
            output_file.write(response.content)
        print(f"File downloaded successfully to {output_path}")
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

#############################

# one option means we'll use the default dates and will use the given workbook name
if len(sys.argv) == 1 or sys.argv[1] == '-h':
    print(f"Pulls PNG images from Tableau Server")
    print(f"Usage:\npython {sys.argv[0]} WORBOOK_NAME [nowcast_week] [weighted_week]")
    exit(0)

if len(sys.argv) > 1:
    workbook_name = sys.argv[1]

#if we have extra options, they should be dates
if len(sys.argv) > 2 and len(sys.argv) <5:
    print("Setting dates")
    query_dates = { 'Nowcast': sys.argv[2] , 'Weighted': sys.argv[3]}
else: 
    query_dates = get_dates()

# Get the token
try:
    auth_token = get_token()
except:
    print("Token not obtained successfully. Please check the script.", file=sys.stderr)
    exit(1)


    # Get the workbook id using the obtained token and workbook name
try:
    workbook_id = get_workbook_id(auth_token, workbook_name)
except:
    print(f"Could not find workbook {workbook_name}!", file=sys.stderr)
    exit(1)

output_dir = 'output_' + str(datetime.now().date())
mkdir(output_dir)


print(f"The id attribute for the workbook with name '{workbook_name}' is: {workbook_id}")
print(f"The Nowcast date is set to {query_dates['Nowcast']}; the Weighted date is set to {query_dates['Weighted']}")
for view in view_names:
    try:
        view_id = get_view_id(auth_token, workbook_id, view)
    except:
        print(f"Could not find view {view}!", file=sys.stderr)
        exit(1)
    for d in query_dates.keys():
        if view == 'Mapforppt':
            regionb = ['USA']
        else:
            regionb = regions
        for region in regionb:
            p = 'vf_USA%20or%20HHSRegion%20Parameter={}&vf_Parameters.Weighted%20or%20Nowcast={}&vf_Selected%20week%20ending={}'.format(region, d, query_dates[d])
            filename = f'{output_dir}/{view}_{d}_for_Region_{region}.png'
            get_image(auth_token, view_id, p, filename)



