def frec(entrada):
  with open("mutinterest210206.freq","a") as salida_fas:
      titulo="seq\t1163\t1167\n"
      salida_fas.write(titulo)
      with open(entrada,"r") as entrada_fastas:
          for line in entrada_fastas:
              line=line.replace("\n","").split("\t")
              sentencia=line[0]+"\t"+line[1][25049-1]+line[1][25050-1]+line[1][25051-1]+"\t"+line[1][25061-1]+line[1][25062-1]+line[1][25063-1]+"\n"
              print(len(line[1]))
              salida_fas.write(sentencia)
              
frec("lin_fa.fasta")
