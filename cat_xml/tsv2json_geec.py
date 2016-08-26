import sys
import csv
import StringIO
import json

TEST = """Name  Age Address
Paul    23  1115 W Franklin
Bessy the Cow   5   Big Farm Way
Zeke    45  W Main St
"""

def make_json(csv_reader):
    rows = [row for row in csv_reader]
    json_content = {"datasets":rows}
    return json.dumps(json_content)

def main():
    if sys.argv[1] == "test":
        csv_file = StringIO.StringIO(TEST)
    else:
        csv_path = sys.argv[1]
        csv_file = open(csv_path, 'r')
    csv_reader = csv.DictReader(csv_file, dialect='excel-tab')
    print make_json(csv_reader)

if __name__ == "__main__":
    main()