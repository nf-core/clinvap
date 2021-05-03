#!/usr/bin/env python3
from docxtpl import DocxTemplate
import sys
import json

json_report = sys.argv[1]
word_template = sys.argv[2]
output_file = sys.argv[3]


doc = DocxTemplate(word_template)

with open(json_report) as f:
    context = json.load(f)

doc.render(context)
doc.save(output_file)