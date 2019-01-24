import sys
sys.path.append("C:/Python27/Lib/site-packages")
import urllib
import getopt
import json
import os
from Bio import SeqIO
from Bio.Blast import NCBIXML

def readFasta(fastaFileName):
    handle = open(fastaFileName)
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return record_dict

def readBlastXML(blastFileName):
    handle = open(blastFileName)
    return NCBIXML.parse(handle)

def getGiNumbers(blastRecords):
    giNumbers=[]
    for b in blastRecords:
        #print(b.query)
        if(not b.alignments):
            #print("No Alignments")
            giNumbers.append('1234567890')
        else:
            giNumbers.append(b.alignments[0].title.split('|')[1])
    return giNumbers

def getGeneSymbols(giNumbers):
    url = "https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input=ginumber&inputValues="
    url += ','.join(giNumbers) + "&outputs=genesymbol"
    u = urllib.urlopen(url)
    response = u.read()
    #f = open("temp_response", 'r')
    #response = f.read()
    #f.write(response)
    data = json.loads(response)
    #pprint(data)
    responseDict={}
    for resp in data:
        responseDict[resp["InputValue"]] = resp["Gene Symbol"]

    return responseDict

def usage():
    print("USAGE: python seqCaptAnnotate.py -b BLAST.XML -f CONSENSUS.fasta -o OUTPUTFILE.fasta")

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hb:f:o:", ["help", "blastFile=", "fastaFile=", "outputFile="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-b", "--blastFile"):
            blastFileName = arg
        elif opt in ("-f", "--fastaFile"):
            fastaFileName = arg
        elif opt in ("-o", "--outputFile"):
            outputFile = arg

    #print(blastFileName, fastaFileName, referenceName)
    fastaDict = readFasta(fastaFileName)
    blastRecords = readBlastXML(blastFileName)
    giNumbers = getGiNumbers(blastRecords)
    blastRecords = readBlastXML(blastFileName)
    geneSymbols = getGeneSymbols(giNumbers)
    
    with open(outputFile, 'w') as f:
        for b in blastRecords:
            if (not b.alignments):#no alignment
                fastaDict[b.query].id = fastaFileName.split('.')[0] \
                    + " No Alignment" + " | Assembled from " \
                    + fastaDict[b.query].id + " | Most Similar to " \
                    + "No Alignment"
            else:#alignment found
                fastaDict[b.query].id = fastaFileName.split('.')[0] \
                    + " " \
                    + geneSymbols[b.alignments[0].title.split('|')[1]].upper()\
                    + " | Assembled from " + fastaDict[b.query].id \
                    + " | Most Similar to " \
                    + geneSymbols[b.alignments[0].title.split('|')[1]].upper()\
                    + b.alignments[0].title
            #For some reason, SeqIO.write() reprints old fastaDict[b.query].id
            #Problem is unresolvable (see end)
            #SeqIO.write(fastaDict[b.query], f, 'fasta')
            f.write(">%s\n" % fastaDict[b.query].id)
            count = 0
            for i in fastaDict[b.query]:
                if (count == 60):
                    f.write("\n")
                    count = 0
                f.write(i)
                count = count + 1
            f.write("\n")
    
    #for b in blastRecords:
    #    if(not b.alignments):
    #        if(len(fastaDict[b.query].id.split('_')) > 1):
    #            fastaDict[b.query].id = fastaFileName.split('.')[0] + " No Alignment" + " | Assembled from " + fastaDict[b.query].id + " | Most Similar to " + "No Alignment"
    #        else:
    #            fastaDict[b.query].id = fastaFileName.split('.')[0] + " No Alignment" + " | Assembled from " + fastaDict[b.query].id + " | Most Similar to " + "No Alignment"
    #    else:
    #        if(len(fastaDict[b.query].id.split('_')) > 1):
    #            fastaDict[b.query].id = fastaFileName.split('.')[0] + " " + geneSymbols[b.alignments[0].title.split('|')[1]].upper() + " | Assembled from " + fastaDict[b.query].id + " | Most Similar to " + geneSymbols[b.alignments[0].title.split('|')[1]].upper() + b.alignments[0].title 
    #        else:
    #            fastaDict[b.query].id = fastaFileName.split('.')[0] + " " +geneSymbols[b.alignments[0].title.split('|')[1]].upper() + " | Assembled from " + fastaDict[b.query].id + " | Most Similar to " + geneSymbols[b.alignments[0].title.split('|')[1]].upper() + b.alignments[0].title 

    #os.remove()
    #f = open(outputFile, 'w')
    #f.close()
    #f = open(outputFile, 'a')
    #for i in fastaDict:
    #    SeqIO.write(fastaDict[i], f, 'fasta')

#             <sequence file name> <top blast hit gene name> | Assembled from
# > <reference sequence name> <reference gene name> | Most Similar to <top
# > blast hit species name> <top blast hit gene name> <top blast hit
# > accession number>

    #print(len(giNumbers))

    #The reason why SeqIO.write() reprints old fastaDict[b.query].id is because
    #of another private member called .description. When fastaDict was created,
    #it used SeqIO.to_dict(), which set .description as well as .id. They both
    #are initialized with the same value (.id = .description). The only
    #difference is that only .id can be modified by the user. When calling
    #SeqIO.write(), other subfunctions from other classes are called, but it
    #eventually reaches FastaWriter.write_record(). This method checks whether
    #.description includes .id at the start of itself. Since .id was modified
    #(see above), this check returns false, leading to a new variable:
    # title = "%s %s" % (id, description). This variable is written using:
    # .write(">%s\n" % title), thus appending description (which is the old .id)
    #, which is unwanted and unavoidable if changing .id.

if __name__ == '__main__':
    main(sys.argv[1:])