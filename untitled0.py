# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 11:42:30 2021

@author: edunu
"""
import requests
from bs4 import BeautifulSoup as bs
import time
import re
import random

def uniprot(name):
    name = name.replace(" ", "+")
    up_query_url = "https://www.uniprot.org/uniprot/?query={}&sort=score"
    up_entry_url = "https://www.uniprot.org/uniprot/{}"
    """
    if not name:
        name = input("Enter protein name: ").replace(" ", "+")
    """
    query_url = up_query_url.format(name)
    seen_hits = 0
    id__ = None
    sequences = []
    while True:
        if id__ is None:
            query_page = requests.get(query_url).text
            time.sleep(3)
            query_up = bs(query_page, "html.parser")
        
            up_hits = query_up.findAll("td", class_="entryID")
            up_total_hits = query_up.find("strong",
                                               class_="queryResultCount")
            print("There are {} hits for this protein name and you have "
                  "seen {}. Please, select UniProt entry:"
                  "".format(up_total_hits.text, seen_hits))
            while True:
                for hit in up_hits:
                    seen_hits += 1
                    row = hit.parent.findAll("td")
                    info = [cell for cell in row]
                    id_ = info[2].text
                    name = info[4].find("div", class_="long").text
                    sp = info[6].text
                    length = info[7].text
                    prompt = ("- Entry {}: {}.\n"
                              "-Species: {}\n-Length: {} aa.\n  ").format(id_, 
                                                                          name, 
                                                                          sp, 
                                                                          length)
                    
                    try:                                                      
                        correct_protein = input(prompt)
                    except KeyboardInterrupt:
                        return sequences
                    
                    if correct_protein != "":
                        id__ = hit.text
                        long_name = name
                        organism = sp
                        hit_url = up_entry_url.format(id__)
                        uniprot_page = requests.get(hit_url).text
                        time.sleep(3)
                        soup = bs(uniprot_page, "html.parser")
                        sequence = soup.find("pre", class_="sequence").text
                        for thing in " 0123456789":
                            sequence = sequence.replace(thing, "")
                        sequence = ">{} | {} | {}\n{}".format(correct_protein,
                                                                   id__,
                                                                   long_name,
                                                                   sequence)
                        sequences.append(sequence)
                    elif seen_hits % 25 == 0:
                        if seen_hits == 25:
                            up_query_url = query_url+"&offset={}"
                        query_url = up_query_url.format(seen_hits)
                        break
                if correct_protein == "q" or seen_hits % 25 == 0: 
                    break

        else:
            break
    return sequences
    
trpa1 = uniprot("trpa1")

chicken = uniprot("trpa1 chicken")
squirrel = uniprot("trpa1 squirrel")

trpa1 = trpa1 + chicken + squirrel

f = open("trpa1_orthologs.fasta", "w")
f.write("\n\n".join(trpa1))
f.close()
#%%
def clustal_align(sequences, email, filename = False):
    """
    Performs an alignment of the required protein sequences using EBI 
    Clustal Omega website. The result is printed to the console and sent to
    the email introduced. User can introduce a list of Protein instances,
    or the user will be prompted to enter the name of the proteins to be
    aligned and indicate the desired UniProt entry.        

    Parameters
    ----------
    sequences : list, optional
        List of Protein instances to be aligned. If False, the user is
        prompted to enter the names of the proteins whose sequences are to
        be aligned. The default is False.
    email : string
        Email where the result of the alignment will be sent. Besides, this
        result will be printed to the console.

    Returns
    -------
    None.

    """

    clustal = "https://www.ebi.ac.uk/Tools/services/web_clustalo/toolform.ebi"

    myobj = {'tool': 'clustalo',
             'stype': "protein",
             'sequence': "\n\n".join(sequences),
             'outfmt': "clustal_num",
             'notification': 'on',
             'email': email,
             'title': "job_" + str(random.randint(0,1e6)),
             'submit': 'Submit'}
    
    session = requests.session()
    x = requests.post(clustal, data = myobj, allow_redirects=True)
    i = 0
    while True:
        try:
            suffix = bs(x.text, "html.parser").find("div", id="ebi_mainContent").a.get("href")
            result_url = "https://www.ebi.ac.uk/Tools/services/web_clustalo/" + suffix
            
            page = requests.get(result_url).text
            soup = bs(page, "html.parser")
            
            alignment = soup.find("pre").text
    
            print(alignment, end = "\n\n")
            
            if filename:
                f = open(filename, "w")
                f.write(alignment)
                f.close()
                print(filename, "has been written.")
                
            break
            
        except:
            print("Nope")
            time.sleep(10)
            if i > 3:
                break
            i += 1
  
with open("Desktop/trpa1/trpa1_orthologs.fasta", "r") as f:
    trpa1 = "".join(f.readlines())

trpa1 = ["\n".join(trpa1.split("\n")[i:i+2]) for i in range(0,36, 3)]
    
clustal_align(trpa1, "edunueeeve@gmail.com", "alignments.txt")
clustal_align(trpa1[0:3] + trpa1[6:], "edunueeeve@gmail.com", "reduced_alignments.txt")



#%%
def get_alignment_dict_and_scores(file):
    f = open(file, "r")
    alignment = f.readlines()
    f.close()
    align = {}
    scores = []
    a = 0
    for line in alignment[1:]:
        try:
            org, seq = re.search("(\w+)\s+(\S+)", line).groups()
            try:
                align[org] += seq
                a = 1
            except KeyError:
                align[org] = seq
                a = 1
        except AttributeError:
            if len(line.strip()) > 10 or a == 1:
                scores.append(line)
                a = 0
        
    score_s = "".join([score[14:74].strip("\n") for score in scores])
    return align, score_s

def compute_positions_and_limits(align, sp, low, high, digits):
    chars = 1
    positions = []
    for i,t in enumerate(list(align[sp])):
        if t in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            
            pos = str(chars)
            while len(pos) < digits:
                pos = " " + pos
            positions.append(pos)
            if chars == low:
                low_limit = i
            elif chars == high:
                high_limit = i
            chars += 1
        else:
            positions.append(" "*digits)
    return positions, low_limit, high_limit

def print_positions(position_list, low, high, margin):
    digits = len(position_list[0])
    position_dict = {}
    for i in range(digits):
        position_dict[str(i)] = [pos[i] for pos in position_list]
    
    last_key = str(digits - 1)
    for n in range(digits):
        print("\n" + " "*margin, end = "")
        for i in range(low, high):
            try:
                if int(position_dict["3"][i]) % 5 == 0:
                    print(position_dict[str(n)][i], end = "")
                else:
                    print(" ", end = "")
            except ValueError:
                print(" ", end = "")
                

def print_alignment(alignment_dict, scores, ref_sp, low, high, line_char = 70):
    pos_dict, low_lim, up_lim = compute_positions_and_limits(alignment_dict, 
                                                             ref_sp, low, high,
                                                             4)
    
    window_len = up_lim - low_lim
    margin = 9
    lines = window_len / line_char
    lines = int(lines + 1) if lines > int(lines) else int(lines)
    beg = low_lim
    end = up_lim if lines == 1 else beg + line_char
    
    for i in range(lines):
        print_positions(pos_dict, beg, end, margin)
        print("")
        
        for k in alignment_dict:
            print(k + " "*(8-len(k)), align[k][beg:end])
        
        print(" "*margin + scores[beg:end])
            
        beg += line_char
        end = end + line_char if i < (lines - 2) else up_lim
    
align, scores = get_alignment_dict_and_scores("alignments.txt")
print_alignment(align, scores, "human", 600, 711, 90)