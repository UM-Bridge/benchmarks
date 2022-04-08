import marko
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename", required=True,
                    help="check benchmark documentation markdown FILE for correctness", metavar="FILE")

args = parser.parse_args()


file = open(args.filename,mode='r')
test_file_content = file.read()
file.close()

test_doc = marko.parse(test_file_content)

undefined_headings = {"Overview": 2, "Purpose": 2, "Run": 2, "Properties": 2, "Configuration": 3, "Description": 3, "Model": 2}

for entry in test_doc.children:
  if entry.get_type() == "Heading":
    heading_text = entry.children[0].children
    if heading_text in undefined_headings:
      print(f"Found heading '{heading_text}' at level {entry.level}; expected at {undefined_headings[heading_text]}")
      assert(entry.level == undefined_headings[heading_text])
      del undefined_headings[heading_text]

print(f"Missing headings: {len(undefined_headings)} ({undefined_headings})")
assert len(undefined_headings) == 0
print("Test passed!")