import random

def coger_random(fichero,numero):
  lista_fichero=[]
  with open(fichero,"r") as entrada:
    for line in entrada:
      line=line.replace("\n","")
      lista_fichero.append(line)
  with open("random"+fichero,"w") as salida:
    randomlist = random.sample(lista_fichero, numero)
    for i in randomlist:
      line=i+"\n"
      salida.write(line)


coger_random("1.txt",6)
coger_random("2.txt",13)
coger_random("3.txt",27)
coger_random("4.txt",1037)
coger_random("5.txt",1001)
coger_random("6.txt",509)
coger_random("7.txt",524)
coger_random("8.txt",442)
coger_random("9.txt",440)
coger_random("10.txt",503)
coger_random("11.txt",889)
coger_random("12.txt",901)
coger_random("13.txt",108)