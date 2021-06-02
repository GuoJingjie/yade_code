#!/usr/bin/env python3

# Create backup of questions/answers of launchpad
# https://askubuntu.com/questions/632964/retrieve-launchpad-answers-programatically
# https://launchpad.net/+apidoc/devel.html#question

from launchpadlib.launchpad import Launchpad
import urllib.request
import json

lp = Launchpad.login_anonymously(
    'Yade launchpad questions  backup', version='devel')
project = lp.projects['yade']

# Get the list of questions
questions = project.searchQuestions()

print(f"We have got {len(questions)} questions")

owner_links = set()
questions_dump = {}

# Go through the list and get data
# Please do not run it too often, otherwise we will be banned by copmany, running LP
for n, i in enumerate(questions):
    print(f"Loading question ({n+1}/{len(questions)}) {i.id}: {i.title}")

    # Get question
    with urllib.request.urlopen(i.self_link) as url:
        data_question = json.loads(url.read().decode())
        owner_links.add(data_question['owner_link'])

    # Get messages in question
    with urllib.request.urlopen(i.messages_collection_link) as url:
        data_messages = json.loads(url.read().decode())
        for d in data_messages["entries"]:
            owner_links.add(d['owner_link'])
    questions_dump[i.id] = {
        'question': data_question, 'messages': data_messages}

# Dump questions
with open(f"questions_dump.json", 'w') as f:
    json.dump(questions_dump, f)


owners = []
for i in owner_links:
    try:
        with urllib.request.urlopen(i) as url:
            data_owner = json.loads(url.read().decode())
            owners.append(data_owner)
    except:
        print(f"Exception occurred with {i}")

# Dump owners into json
with open(f"owners.json", 'w') as f:
    json.dump(owners, f)
