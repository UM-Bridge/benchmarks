import glob
import marko
import re
from argparse import ArgumentParser

for filename in glob.glob('../benchmarks/*/README.md') + glob.glob('../models/*/README.md'):
  print(f"Checking file {filename}")

  file = open(filename,mode='r')
  test_file_content = file.read()
  file.close()

  test_doc = marko.parse(test_file_content)

  undefined_headings = {"Overview": 2, "Authors": 2, "Run": 2, "Properties": 2, "Mount directories": 2, "Source code": 2, "Description": 2}

  for entry in test_doc.children:
    if entry.get_type() == "Heading":
      heading_text = entry.children[0].children
      if heading_text in undefined_headings:
        print(f"Found heading '{heading_text}' at level {entry.level}; expected at {undefined_headings[heading_text]}")
        assert(entry.level == undefined_headings[heading_text])
        del undefined_headings[heading_text]

  print(f"Missing headings: {len(undefined_headings)} ({undefined_headings})")

  h1_count = 0
  errors = []
  lines = test_file_content.splitlines()
  for i, line in enumerate(lines):
      if re.match(r'^# ', line):
          h1_count += 1
          if h1_count > 1:
              errors.append(f"Line {i+1}: More than one H1 header found: {line.strip()}")

  if errors:
      print("\n".join(errors))

  assert len(undefined_headings) == 0
  print("Test passed!")
  print("")
