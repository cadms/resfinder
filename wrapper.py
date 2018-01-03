#!/usr/local/anaconda/bin/python2.7
# ########################################################################### #
# #########                    CGE Service Wrapper                    ####### #
# ########################################################################### #
# This script is part of the CGE Service structure
# --> The input/arguments are validated
# --> Create service folder and subfolders
# --> The uploaded files are copied to the 'Upload' directory
# --> Log service submission in SQL
# --> Setup and execution of service specific programs
# --> The success of the service is validated
# --> A HTML output page is created if in webmode
# --> Downloadable files are copied/moved to the 'Download' Directory
# --> The SQL service entry is updated
# --> The files are zipped for long-term storage
import sys
import os
import time
import random
import subprocess
import re
from string import Template

# INCLUDING THE CGE MODULES (No need to change this part)
sys.path.append("/home/data1/services/CGEpipeline/CGEmodules")
from assembly_module_2 import (PrepareAssembly, MakeAssembly,
                               printContigsInfo, AddAssemblyToProgList)
from functions_module_2 import (printDebug, copyFile, program,
                                createServiceDirs, getArguments, paths,
                                makeFileList, fileUnzipper,
                                printInputFilesHtml, fileZipper,
                                PrepareDownload, UpdatePaths, proglist,
                                CheckFileType, printOut, moveFile, setDebug,
                                add2DatabaseLog, dlButton, GraceFullExit,
                                FixPermissions, tsv2html)

#  Modules used to create the extended ResFinder output (phenotype output)
from phenotype2genotype.isolate import Isolate
from phenotype2genotype.res_profile import PhenoDB
from phenotype2genotype.res_sumtable import ResSumTable

# ########################################################################### #
# #########                         FUNCTIONS                       ######### #
# ########################################################################### #


def create_tab_acquired(isolate, phenodb):
   """ Alternative method to create the downloadeable tabbed result file. This
       method will include the additional information from the phenotype
       database.
   """
   output_str = ("Resistance gene\tIdentity\tAlignment Length/Gene Length\t"
                 "Position in reference\tContig\tPosition in contig\tPhenotype"
                 "\tClass\tPMID\tAccession no.\tNotes\n")

   for unique_id in isolate:
      for feature in isolate[unique_id]:

         # Extract phenotypes
         phenotype_out_list = []
         phenotype = phenodb[feature.unique_id]

         # Append stars to phenotypes that are suggested by the curators and
         # not published
         for antibiotic in phenotype.phenotype:
            if(antibiotic in phenotype.sug_phenotype):
               antibiotic = antibiotic + "*"
            phenotype_out_list.append(antibiotic)

         phenotype_out_str = ",".join(phenotype_out_list)

         output_str += (feature.hit.name + "\t"
                        + str(feature.hit.identity) + "\t"
                        + str(feature.hit.match_length)
                        + "/" + str(feature.hit.ref_length) + "\t"
                        + str(feature.hit.start_ref)
                        + ".." + str(feature.hit.end_ref) + "\t"
                        + feature.seq_region + "\t"
                        + str(feature.start)
                        + ".." + str(feature.end) + "\t"
                        + phenotype_out_str + "\t"
                        + ",".join(phenotype.ab_class) + "\t"
                        + ",".join(phenotype.pmid) + "\t"
                        + feature.hit.acc + "\t"
                        + phenotype.notes + "\n")

   # Find AMR classes with no hits
   no_class_hits = []
   for ab_class in phenodb.antibiotics:
      if(ab_class not in isolate.resprofile.resistance_classes):
         no_class_hits.append(ab_class)

   if(no_class_hits):
      output_str += ("\nNo hits found in the classes: "
                     + ",".join(no_class_hits) + "\n")

   return output_str


def create_html_acquired(isolate, phenodb, indent="  "):
   """ Alternative method to the method tab2html_acquired for table creation.
       The idea is that the isolate object already contains all the traditional
       ResFinder information but also the 'new' phenotype information.
   """

   match_cat_translation = {1: "bggrey",
                            2: "bglightgreen",
                            3: "bggreen"}

   class_table_beg = Template(
       "$indent<div class='resresults'>\n"
       "$indent<table width='100%'\n"
       "$indent  <thead><tr><th colspan='10'>$ab_class</th></tr></thead>\n"
       "$indent  <tr>\n"
       "$indent    <th>Resistance gene</th>\n"
       "$indent    <th>Identity</th>\n"
       "$indent    <th>Alignment Length/Gene Length</th>\n"
       "$indent    <th>Position in reference</th>\n"
       "$indent    <th>Contig</th>\n"
       "$indent    <th>Position in contig</th>\n"
       "$indent    <th>Phenotype</th>\n"
       "$indent    <th>PMID</th>\n"
       "$indent    <th>Accession no.</th>\n"
       "$indent    <th>Notes</th>\n"
       "$indent  </tr>\n"
   )

   class_table_row = Template(
       "$indent  <tr class='$match_class'>\n"
       "$indent    <td style='text-align:center'>$gene</td>\n"
       "$indent    <td style='text-align:center'>$identity</td>\n"
       "$indent    <td style='text-align:center'>$al_length/$ref_length</td>\n"
       "$indent    <td style='text-align:center'>$start_ref..$end_ref</td>\n"
       "$indent    <td style='text-align:center'>$seq_region</td>\n"
       "$indent    <td style='text-align:center'>$start..$end</td>\n"
       "$indent    <td style='text-align:center'>$pheno</td>\n"
       "$indent    <td style='text-align:center'>$pmid</td>\n"
       "$indent    <td style='text-align:center'>"
       "<a href='http://www.ncbi.nlm.nih.gov/nuccore/$acc'> $acc </a></td>\n"
       "$indent    <td style='text-align:center'>$note</td>\n"
       "$indent  </tr>\n"
   )

   output_html = ""

   for ab_class in isolate.resprofile.resistance_classes:
      output_html += class_table_beg.substitute(ab_class=ab_class,
                                                indent=indent)

      for feature in isolate.resprofile.resistance_classes[ab_class]:

         # Extract phenotypes
         phenotype_out_list = []
         phenotype = phenodb[feature.unique_id]

         # Append stars to phenotypes that are suggested by the curators and
         # not published
         for antibiotic in phenotype.phenotype:
            if(antibiotic in phenotype.sug_phenotype):
               antibiotic = antibiotic + "*"
            phenotype_out_list.append(antibiotic)

         phenotype_out_str = ", ".join(phenotype_out_list)

         output_html += class_table_row.substitute(
             indent=indent,
             match_class=match_cat_translation[feature.hit.match_category],
             gene=feature.hit.name,
             identity=feature.hit.identity,
             al_length=feature.hit.match_length,
             ref_length=feature.hit.ref_length,
             start_ref=feature.hit.start_ref,
             end_ref=feature.hit.end_ref,
             seq_region=feature.seq_region,
             start=feature.start,
             end=feature.end,
             pheno=phenotype_out_str,
             pmid=", ".join(phenotype.pmid),
             acc=feature.hit.acc,
             note=phenotype.notes)

      output_html += indent + "</table>\n"
      output_html += indent + "</div>\n"
      output_html += indent + "<br><br>\n\n"

   # Find AMR classes with no hits
   no_class_hits = []
   for ab_class in phenodb.antibiotics:
      if(ab_class not in isolate.resprofile.resistance_classes):
         no_class_hits.append(ab_class)

   if(len(no_class_hits) > 0):
      output_html += (indent + "<strong><p>No hits were found for the "
                      "following classes:</p></strong><br>\n")
      output_html += indent + "<ul>\n"
      for ab_class in no_class_hits:
         output_html += indent + "  <li>" + ab_class + "</li>\n"
      output_html += indent + "</ul>\n"
      output_html += indent + "<br><br>\n\n"

   return output_html


def tab2html_acquired(fp):
   printDebug('Writing HTML table for %s' % fp)

   if os.path.exists(fp):
      with open(fp, 'r') as f:
         intable = False
         tsize = 0
         ttitle = ''

         for i in f:
            i = i.strip()

            if i == '':
               # Ending the table and making a new line
               printOut("      </table>")
               printOut("<br>")
               intable = False
            elif i.startswith('Resistance') and intable:
               tmp = i[0:].split('\t')
               tsize = len(tmp)
               # Print Header Line
               printOut("        <tr><th>%s</th></tr>" %
                        ("</th><th>".join(tmp)))
            elif intable:
               td = i.split("\t")
               # If the line are not tab seperated the line is just printed
               if (len(td) < 2):
                  printOut("        <tr class='bgred'><td "
                           "style='text-align:center'>%s</td></tr>" % (i))
               else:
                  ID = td[1]
                  tmp = td[2].split('/')
                  HSP = int(float(tmp[1]))
                  query = int(float(tmp[0]))
                  acc = td[-1]
                  td[-1] = ("<a href='http://www.ncbi.nlm.nih.gov/nuccore/%s'>"
                            " %s </a>" % (acc, acc))

                  # Setting the background color
                  if (ID == '100.00') and (HSP == query):
                     printOut("<tr class='bggreen'>")
                  elif (ID != '100.00' and (HSP == query)):
                     printOut("<tr class='bglightgreen'>")
                  elif (ID == '100.00' and (HSP != query)):
                     printOut("        <tr class='bggrey'>")
                  else:
                     printOut("        <tr class='bgred'>")

                  # Printing the result
                  printOut("          <td style='text-align:center'>%s</td>" %
                           ("</td><td style='text-align:center'>".join(td)))
            else:
               # Start table
               intable = True

               # Setting the table class
               printOut("      <table class='center virresults' width='100%'>")
               printOut("        <tr><th colspan='8'>%s</th></tr>" % (i))

         # END TABLE IF OPEN
         if intable:
            printOut("</table>")


def printExtendOutput_acquired(query, subject):
   #Define class
   #print("<style type='text/css'> .bglightgreen{background-color:#BEDCBE;} .bggrey{background-color:#9D9D9D;}.bggreen{background-color:#6EBE50;}.bgred{background-color:#BF2026;}.virresults td, th { border: 2px solid #FFFFFF; padding: 4px 16px 1px; text-align: left; vertical-align: middle;} .virresults th,td{font-size: 14px;padding: 4px 5px 1px 5px;text-align: center;} .virresults th{color:#FFF; background-color:#3956A6;}.virresults td a{color:#3956A6;}</style>")
   printOut("<style type='text/css'> .bggrey{background-color: #AAAAAA;} .bgred{background-color: #BF2026;} .bggreen{background-color: #6EBE50;} tr.header{font-weight: bold; text-align: center; background-color: grey;} .w70{width:70%;} span, p{font-size: 16px;}</style>")


   # Print extended output button
   printOut("<script type='text/javascript'>function printOutput(tag){var ele = document.getElementById(tag).style;if (ele.display=='block'){ele.display='none';}else{ele.display='block';}}</script>\n")
   printOut("<center><button type='button' onclick='printOutput(&quot;eo&quot;)'>extended output</button></center><div id='eo' class='hide'>")

   # Start output
   printOut("<pre><br>")

   # Get output sequences
   if os.path.exists(query) and os.path.exists(subject):
      with open(query, 'r') as q, open(subject, 'r') as s:
         query_sequences = dict()
         subject_sequences = dict()
         nrLine = 0

         # Reading in the reference lines
         slines = s.readlines()
         q_header = ''
         for ql in q:
            ql = ql.strip("\n")
            sl = slines[nrLine].strip("\n")
            if ql.startswith(">"):
               if q_header:
                  if s_header in subject_sequences:
                     query_sequences[s_header].append((q_header, q_seq))
                  else:
                     query_sequences[s_header] = [(q_header, q_seq)]
                     subject_sequences[s_header] = s_seq
               # Saving new header and intilising seqeunce
               q_header = ql[1:]
               s_header = sl[1:]
               q_seq = ''
               s_seq = ''
            else:
               q_seq = q_seq + ql
               s_seq = s_seq + sl
            nrLine += 1
         # Saving the last sequence
         if s_header in subject_sequences:
            query_sequences[s_header].append((q_header, q_seq))
         else:
            query_sequences[s_header] = [(q_header, q_seq)]
            subject_sequences[s_header] = s_seq

      # Predefining variables
      startPos = 0
      endPos = 0
      curPos = 0
      ID = 0
      HSP = 0
      ref = 0

      # Looping over saved hits
      for s_header in subject_sequences:
         if len(query_sequences[s_header]) > 1:
            ID = list()
            HSP = list()
            ref = list()
            q_lines = list()
            s_lines = list()
            q_headers = list()
            curPos = list()

            for hit in xrange(len(query_sequences[s_header])):
               q_header = query_sequences[s_header][hit][0]
               q_headers.append(q_header)
               match = re.search('.*ID:\s*(\d+\.*\d*).*Alignment Length/Gene Length:\s+(\d+)/(\d+).*', q_header)
               ID.append(match.group(1))
               HSP.append(int(match.group(2)))
               ref.append(int(match.group(3)))
               q_lines.append(query_sequences[s_header][hit][1])
               s_lines.append(subject_sequences[s_header])
               curPos.append(0)
               printOut(("<span style='font-size:15px; color:grey'>Hit %s: %s</span><br>")%((hit + 1), q_header))

            # Printing the output depending on it being a perfect match or not
            s_lines = subject_sequences[s_header]
            for i in range(0, len(s_lines), 60):
               s_chars = s_lines[i:i + 60]
               finalRefLine = (("<span style='font-size:14px; color:black'>Resistance gene seq:   </span><span>%s</span>")%(s_chars))
               finalHitLine = ''

               for hit in xrange(len(query_sequences[s_header])):
                  q_chars = q_lines[hit][i:i + 60]
                  if (ID[hit] == '100.00') and (HSP[hit] == ref[hit]):
                     printOut(("<span style='font-size:14px; color:black'>Hit in genome %s:       </span><span class='bggreen'>%s</span><br><br>")%((hit + 1),q_chars))
                  else:
                     m = re.search('.*Positions\sin\sreference:\s*(\d+)..(\d+).*', q_headers[hit])
                     startPos = int(float(m.group(1)))
                     endPos = int(float(m.group(2)))
                     finalHitLine += "<br><span style='font-size:14px; color:black'>Hit in genome %s:       </span>"%((hit + 1))

                     # Going through every base and defining the color
                     pos = 0
                     for char in q_chars:
                        skip = 0
                        # Getting the reference base unless the query is longer than the reference
                        if pos > (len(s_chars) - 1):
                           pos += 1
                           skip = 1
                        else:
                           sChar = s_chars[pos]

                        # If the posistion is not included in the alignment the background is gray
                        # if the query is longer than the refernece then the posisions are not inlcuded
                        if skip == 1:
                           pass
                        elif (curPos[hit] < (startPos - 1)) or (curPos[hit] > (endPos - 1)):
                           color = 'grey'
                        else:
                           # If the bases align the color is green ortherwise it is red
                           if char == sChar:
                                 color = 'green'
                           else:
                                 color = 'red'

                        # Adding the base with the right color
                        finalHitLine += "<span class='bg%s'>%s</span>"%(color, char)
                        pos += 1
                        curPos[hit] += 1

               # Printing the alignmet lines.
               printOut(('%s%s<br>')%(finalRefLine, finalHitLine))
            printOut("<br><br>")
         else:
            q_header = query_sequences[s_header][0][0]
            q_line = query_sequences[s_header][0][1]
            s_line = subject_sequences[s_header]

            # Getting inforamiton from query header
            match = re.search('.*ID:\s*(\d+\.*\d*).*Alignment Length/Gene Length:\s+(\d+)/(\d+).*', q_header)
            ID = match.group(1)
            HSP = int(float(match.group(2)))
            ref = int(float(match.group(3)))

            # Printing the output depending on it being a perfect match or not
            if (ID == '100.00') and (HSP == ref):
               printOut(("<span style='font-size:15px; color:grey'>%s</span><br>")%(q_header))
               for i in range(0, len(q_line), 60):
                  q_chars = q_line[i:i + 60]
                  s_chars = s_line[i:i + 60]
                  printOut(("<span style='font-size:14px; color:black'>Resistance gene seq: </span><span class='bggreen'>%s</span>")%(q_chars))
                  printOut(("<span style='font-size:14px; color:black'>Hit in genome:       </span><span class='bggreen'>%s</span><br><br>")%(s_chars))
            else:
               printOut(("<span style='font-size:15px; color:grey'>%s</span><br>")%(q_header))
               m = re.search('.*Positions\sin\sreference:\s*(\d+)..(\d+).*', q_header)
               startPos = int(float(m.group(1)))
               endPos = int(float(m.group(2)))
               curPos = 0

               for i in range(0, len(q_line), 60):
                  q_chars = q_line[i:i + 60]
                  s_chars = s_line[i:i + 60]
                  # Creating a string with the html code
                  finalVirLine = "<br><br><span style='font-size:14px; color:black'>Resistance gene seq: </span>"
                  finalHitLine = "<br><span style='font-size:14px; color:black'>Hit in genome:       </span>"

                  # Going through every base and defining the color
                  pos = 0
                  for char in q_chars:
                     skip = 0
                     # Getting the reference base unless the query is longer than the reference
                     if pos > (len(s_chars) - 1):
                        pos += 1
                        skip = 1
                     else:
                        sChar = s_chars[pos]

                     # If the posistion is not included in the alignment the background is gray
                     # if the query is longer than the refernece then the posisions are not inlcuded
                     if skip == 1:
                        pass
                     elif (curPos < (startPos - 1)) or (curPos > (endPos - 1)):
                        color = 'grey'
                     else:
                        # If the bases align the color is green ortherwise it is red
                        if char == sChar:
                              color = 'green'
                        else:
                              color = 'red'

                     # Adding the base with the right color
                     finalVirLine += "<span class='bg%s'>%s</span>"%(color, sChar)
                     finalHitLine += "<span class='bg%s'>%s</span>"%(color, char)
                     pos += 1
                     curPos += 1

                  # Printing the alignmet lines.
                  printOut(('%s%s')%(finalVirLine,finalHitLine))
         printOut("<br><br>")
      # End output
      printOut("</pre><br><br></div>")

def tab2html_point(tab):
  printOut("<style type='text/css'>.bggreen{background-color:#6EBE50;} /*PERFECT HIT COLOR*/.bgred{background-color:#BF2026;} /*NO HITS COLOR*/.results th,td{ /*PROPERTIES OF TABLE CELLS*/font-size: 14px;padding: 4px 5px 1px 5px;text-align: center;}.results th{ /*HEADER COLORS*/color:#FFF;background-color:#4169E1;}.results td a{color:#4169E1;} /*LINK COLOR*/</style>")

  printDebug('Writing HTML table for %s'%tab)
  if os.path.exists(tab):
     with open(tab, 'r') as f:
        intable = False
        for l in f:
           l = l.strip()
           if re.search('^Chromosomale point', l) != None:
              printOut("<h2> %s </h2>"%(l))
           elif re.search('^Species', l) != None:
              td = l.split(":")
              printOut("<h3>%s:<i class='grey'>%s</i></h3>"%(td[0], td[1]))
           elif l == '':
              # Ending the tabe and making a new line
              if intable: printOut("</table>\n")
              printOut('<br>')
              intable = False
           elif re.search('^Mutation', l) != None and intable:
              tmp = l[0:].split('\t')
              # Print Header Line
              printOut("<tr class='bggrey'><th>%s</th></tr>\n"%("</th><th>".join(tmp)))
           elif intable:
              td = l.split('\t')
              # It the line are not tab seperated the line is just printed
              if (len(td) < 2):
                 printOut("<tr class='bggrey'><td style='text-align:center'>%s</td></tr>\n"%(l))
              else:
                 acc = td[-1]
                 td[-1] = "<a href='http://www.ncbi.nlm.nih.gov/pubmed/%s'> %s </a>"%(acc, acc)
                 printOut("<tr class='bggreen'><td style='text-align:center'>%s</td></tr>\n"%("</td><td style='text-align:center'>".join(td)))
           else:
              # Start table
              intable = True
              # Setting the table class
              printOut("<table class='center results' width='100%'>\n\n")
              printOut("<tr><th colspan='5'>%s</th></tr>\n"%(l))
        # END TABLE IF OPEN
        if intable: printOut("</table>\n")


# ########################################################################### #
# #########                           MAIN                            ####### #
# ########################################################################### #

# SET GLOBAL VARIABLES
setDebug(False)
service, version = "ResFinder", "4.0"

# PARSE ARGUMENTS
args = getArguments('''
selectionbox    technology   -t      VALUE
selectionbox    threshold    -T      VALUE
selectionbox    anti         -a      VALUE
selectionbox    minlength    -L      VALUE
selectionbox    species      -S      VALUE
checkbox        point        -c      VALUE
checkbox        acquired     -acq    VALUE
selectionsbox   showmut      -U      VALUE
''', allowcmd=True)


#----DELETE THIS----#
# VALIDATE REQUIRED ARGUMENTS
if args.technology == None: GraceFullExit("Error: No technology platform was chosen!\n")
if args.uploadPath is None or len(args.uploadPath) < 5 or not os.path.exists(args.uploadPath):
   GraceFullExit("Error: No valid upload path was provided! (%s)\n"%(args.uploadPath))
elif args.uploadPath[-1] != '/': args.uploadPath += '/' # Add endslash
#-------------------#

# SET RUN DIRECTORY (No need to change this part)
if args.reuseid:
  #/srv/www/secure-upload/isolates/5_16_4_2015_162_299_423438//0/
  runID = [x for x in args.uploadPath.split('/') if x!=''][-2]
else:
  runID = time.strftime('%w_%d_%m_%Y_%H%M%S_')+'%06d'%random.randint(0,999999)
paths.serviceRoot = '{programRoot}IO/%s/'%(runID)
paths.isolateRoot = '{programRoot}'
paths.add_path('uploads', '{serviceRoot}Uploads/')

# SET AND UPDATE PATHS (No need to change this part)
UpdatePaths(service, version, '', '', args.webmode)
databases = "/home/data1/services/%s/database_repository"%(service)
blast = "/home/data1/tools/bin/blastn"

# CREATE SERVICE DIRECTORIES (No need to change this part)
createServiceDirs()
stat = paths.Create('uploads')

# LOG SERVICE SUBMISSION IN SQL (No need to change this part)
add2DatabaseLog(service+'-'+version, runID, args.usr, args.ip, args.technology)

# MOVE UPLOADS FROM APP- TO ISOLATE UPLOAD DIRECTORY (No need to change this part)
if stat:
   # Move the uploaded files to the upload directory
   moveFile(args.uploadPath+'*', paths['uploads'])
   # Unzipping uploaded files if zipped
   fileUnzipper(paths['uploads'])
   # GET INPUT FILES from input path
   inputFiles = makeFileList(paths['uploads'])
else:
   GraceFullExit("Error: Could not create upload directory!\n")

#
#----CHANGE THIS----#
#
# GET CONTIGS (No need to change this part UNLESS you dont need assembly)
if args.technology != 'Assembled_Genome':
   # ADD THE ASSEMBLER to the program list
   AddAssemblyToProgList()
   # Prepare the Assembly program for execution
   PrepareAssembly(args.technology, inputFiles)
   # Assemble the reads into contigs
   n50 = MakeAssembly(args.technology)
   # The success of the assembler is validated
   status = progs['Assembler'].GetStatus()
   if status != 'Done' or not isinstance(n50, int):
      GraceFullExit("Error: Assembly of the inputfile(s) failed!\n"%(len(inputFiles)))
else:
   # Validate that only 1 file was submitted
   if len(inputFiles) != 1:
      GraceFullExit("Error: Invalid number of contig files (%s)\n"%(len(inputFiles)))
   # Validate that the uploaded file is fasta
   if CheckFileType(inputFiles[0]) != 'fasta':
      GraceFullExit("Error: Invalid contigs format (%s)!\nOnly the fasta format is recognised as a proper contig format.\n"%(CheckFileType(inputFiles[0])))
   # Add contigsPath
   paths.add_path('contigs', paths['inputs']+'contigs.fsa')
   # Copy file to Input directory
   copyFile(inputFiles[0], paths['contigs'])
#
#-------------------#
#

# ADDING PROGRAMS
if (args.acquired):
  name = "RF"
  workDir = "%s/%s"%(paths.serviceRoot,name)
  os.mkdir(workDir)
  resfinder = program(
   name=name, path=paths['scripts']+'ResFinder-3.0.py', timer=0,
   ptype='python', toQueue=True, wait=False, workDir=workDir,
   args=['-i', paths['contigs'],
         '-t', args.threshold,
         '-l', args.minlength,
         '-d', args.anti,
         '-o', paths['outputs'],
         '-p', databases,
         '-b', blast])
  resfinder.Execute(True)
  proglist.Add2List(resfinder)

if (args.point):
  name = "PF"
  workDir = "%s/%s"%(paths.serviceRoot,name)
  os.mkdir(workDir)
  pointmutations = program(
    name=name, path=paths['scripts']+'PointFinder-3.0.pl', timer=0,
    ptype='perl', toQueue=True, wait=False, workDir=workDir,
    args=['-i', paths['contigs'],
          '-k', args.threshold,
          '-l', args.minlength,
          '-d', args.anti,
          '-s', args.species,
          '-u', args.showmut,
          '-R', paths.serviceRoot])
  pointmutations.Execute(True)
  proglist.Add2List(pointmutations)

# EXECUTION OF THE PROGRAMS
for progname in proglist.list:
  if proglist[progname].status == 'Executing':
     proglist[progname].WaitOn(silent=False)

for progname in proglist.list:
  status = proglist[progname].GetStatus()
  if status != 'Done': GraceFullExit("Error: Execution of the program failed!\n")


# # # CREATE SUMMARY TABLE # # #

pretty_names = {
   "campylobacter": "Campylobacter spp",
   "c.jejuni": "C. jejuni",
   "c.coli": "C. coli",
   "e.faecalis": "E. faecalis",
   "e.faecium": "E. faecium",
   "e.coli": "E. coli",
   "salmonella": "Salmonella spp",
   "s.aureus": "S. aureus"
}
species_name = pretty_names.get(args.species, args.species)

if(args.acquired):
   # Load genotype to phenotype database
   input_for_pheno_db = ("/home/data1/services/%s/database_pheno/"
                         "interpreter_db.txt" % (service))
   res_pheno_db = PhenoDB(input_for_pheno_db)

   # Isolate object stores results
   isolate = Isolate(name=inputFiles[0])

   # isolate.load_resfinder_tab(paths['outputs'] + 'results_table.txt')
   isolate.load_resfinder_tab(paths['outputs'] + 'results_tab.tx')
   isolate.calc_res_profile(res_pheno_db)

   # Create and write the downloadable tab file
   pheno_profile_str = isolate.profile_to_str_table(with_header=True)

   with open(paths['outputs'] + 'pheno_table.txt', 'w') as fh:
      fh.write(pheno_profile_str)

   # Load AMR panels
   input_amr_panels = ("/home/data1/services/%s/database_pheno/"
                       "amr_panels.txt" % (service))
   res_sum_table = ResSumTable(pheno_profile_str)
   res_sum_table.load_amr_panels(input_amr_panels)

   #
   # Create ResFinder Summary table HTML output
   #

   # Entire table is wrapped in the class res-sum and id res-sum-page-wrap
   html_sum_table = '   <div class="res-sum" id="res-sum-page-wrap">\n\n'

   # Extract all panel names and create and id for each to be used in the
   # webpage code
   html_res_sum_panels = []
   html_res_sum_panel_ids = []

   if(species_name in res_sum_table.panels):
       html_res_sum_panels.append(species_name)
       html_res_sum_panel_ids.append("res-sum-0")

   # Create tabs/buttons
   html_sum_table += '     <div class="tab">\n'

   res_sum_tab = Template('      '
                          '<button class="tablinks" '
                          'onclick="openPanel(event, \'$panel_name\')" '
                          'id="defaultOpen">$panel_name</button>\n')

   for panel_name in html_res_sum_panels:
      html_sum_table += res_sum_tab.substitute(panel_name=panel_name)

   html_sum_table += res_sum_tab.substitute(panel_name='Complete')
   html_sum_table += '     </div>\n\n'

   # Create a table for each panel
   if(species_name in res_sum_table.panels):
      html_sum_table += res_sum_table.get_html_panel_table(species_name, 'res-sum-0') + '\n\n'

   html_sum_table += res_sum_table.get_html_panel_table(
      'Complete', 'res-sum-complete') + '\n\n'

   html_sum_table += '   </div>\n\n'

   # WRITE TAB SEPERATED RESFINDER RESULTS
   acquired_tab_str = create_tab_acquired(isolate, res_pheno_db)
   with open(paths['outputs'] + 'resfinder_tab.txt', "w") as fh:
      fh.write(acquired_tab_str)

# # # PRINT HTML # # #

printOut('<h1>%s-%s Server - Results</h1>' % (service, version))

# PRINT THE PROGRAM PARAMETERS TO THE HTML OUTPUT FILE
if len(inputFiles) > 0:
   printInputFilesHtml(inputFiles)

# PRINT THE CONTIGS INFORMATION TO THE HTML OUTPUT FILE
if args.technology != 'Assembled_Genome':
   printContigsInfo(runID=runID)

# PRINT RESFINDER SUMMARY TABLE AND BUTTONS
if(args.acquired):
   printOut(html_sum_table)

   # DOWNLOAD BUTTONS FOR SUMMARY TABLE
   PrepareDownload(paths['outputs'] + 'pheno_table.txt')
   dlButton('Download phenotype table (txt)', 'pheno_table.txt')

   printOut('<br><br>')
   printOut('<h2>Acquired antimicrobial resistance gene - Results</h2>')
   printOut('<br><br>')
   html_acquired_table = create_html_acquired(isolate, res_pheno_db)
   printOut(html_acquired_table)

   PrepareDownload(paths['outputs'] + 'resfinder_tab.txt')
   dlButton('Results as tabseperated file', 'resfinder_tab.txt')

# SCRIPT FOR SUMMARY TABLE
printOut("    <script type=\"text/javascript\" src=\"https://cge.cbs.dtu.dk/services/ResFinder-3.0_Test/cgi-bin/res-sum.js\"></script>")
printOut("    <script>")
printOut("      document.getElementById(\"defaultOpen\").click();")
printOut("    </script>")

printOut("    <br><br>\n")

if (args.acquired):
  # Print output button
  printOut("<script type='text/javascript'> "+
           "function printFunction(id) {"+
           "var e = document.getElementById(id); "+
           "if(e.style.display == 'inline') "+
           "e.style.display = 'none'; "+
           "else "+
           "e.style.display = 'inline'; "+
           "}</script>")
  printOut("<center><button type='button' onclick='printFunction(&quot;eo1&quot;)'> "+
           "Show Acquired antimicrobial resistance results</button></center> "+
           "<div id='eo1' class='hide'>")

  # Start output
  printOut('<h2>Acquired antimicrobial resistance gene - Results</h2>')
  # PREPARE THE DOWNLOADABLE FILE(S)
  _ = PrepareDownload(paths['serviceRoot']+'outputs/results.txt')
  _ = PrepareDownload(paths['serviceRoot']+'outputs/Hit_in_genome_seq.fsa')
  _ = PrepareDownload(paths['serviceRoot']+'outputs/Resistance_gene_seq.fsa')
  _ = PrepareDownload(paths['serviceRoot']+'outputs/results_tab.txt')

  # PRINT TAB FILE(S) AS RESULT TABLES TO THE HTML OUTPUT FILE
  tab=paths['outputs']+'results_table.txt'
  if os.path.exists(tab):
    tab2html_acquired(tab)
    printOut('<br>\n')

  # GETTING THE MINIMUM LENGTH FOR OUTPUT PRINT
  minlength = float(args.minlength) * 100
  threshold = float(args.threshold) * 100

  # PRINT THE SETTINGS
  printOut("<h2>Selected %%ID threshold: &nbsp;<i class='grey'>%d %%</i></h2>"%(threshold))
  printOut("<h2>Selected minimum length: &nbsp;<i class='grey'>%d %%</i></h2>"%(minlength))

  ## PRINT FILE DOWNLOAD-BUTTON(S) TO THE HTML OUTPUT FILE
  dlButton('Results as text', 'results.txt')
  dlButton('Hit in genome sequences', 'Hit_in_genome_seq.fsa')
  dlButton('Resistance gene sequences', 'Resistance_gene_seq.fsa')
  dlButton('Results as tabseperated file', 'results_tab.txt')

  ## PRINT THE EXTENDED OUTPUT
  ref=paths['downloads']+'Resistance_gene_seq.fsa'
  hit=paths['downloads']+'Hit_in_genome_seq.fsa'
  if os.path.exists(ref) and os.path.exists(hit):
    printExtendOutput_acquired(hit, ref)
    printOut('<br>\n')
  # End output
  printOut("<br><br></div>")

if (args.point):
  # Print output button
  printOut("<script type='text/javascript'> "+
           "function printFunction(id) { "+
           "var e = document.getElementById(id); "+
           "if(e.style.display == 'inline') "+
           "e.style.display = 'none'; "+
           "else "+
           "e.style.display = 'inline'; "+
           "}</script>")
  printOut("<center><button type='button' onclick='printFunction(&quot;eo2&quot;)'> "+
           "Show Point mutation results</button></center> "+
           "<div id='eo2' class='hide'>")

  # Start output
  printOut("<br>")
  ## PRINT TAB FILE(S) AS RESULT TABLES TO THE HTML OUTPUT FILE
  point=paths['downloads']+'point_results_table.txt'
  if os.path.exists(point):
    tab2html_point(point)
    printOut('<br>\n')

  ## PRINT FILE DOWNLOAD-BUTTON(S) TO THE HTML OUTPUT FILE
  _ = PrepareDownload(paths['serviceRoot']+'downloads/point_results_tab.txt')
  dlButton('Results as tabseperated file', 'point_results_tab.txt')
  printOut("<br><br></div>")

################################################################################
# LOG THE TIMERS (No need to change this part)
proglist.PrintTimers()

# FIX THE PERMISSIONS OF THE SERVICE ROOT
FixPermissions()

# INDICATE THAT THE WRAPPER HAS FINISHED (No need to change this part)
printDebug("Done")

# GZIP ALL FILES IN THE SERVICE DIRECTORY AND SUBDIRS FOR LONG-TERM STORAGE (No need to change this part)
fileZipper(paths['serviceRoot'])
