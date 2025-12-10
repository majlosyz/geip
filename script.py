from Bio import SeqIO
rekord = SeqIO.read("/home/mkrysinski/results_bakta/bakta.dc/assembly.gbff", "genbank")
liczba_cds = 0
liczba_funkcja_znana = 0
liczba_funkcja_nieznana = 0
dlugosci = []

for cecha in rekord.features:
    if cecha.type == "CDS":
        liczba_cds += 1
        produkt = cecha.qualifiers.get("product", [""])[0].lower()
        
        if produkt == "" or "hypothetical protein" in produkt:
            liczba_funkcja_nieznana += 1
        else:
            liczba_funkcja_znana += 1
        
        dlugosc = len(cecha.location)
        dlugosci.append(dlugosc)

dlugosci.sort()
n = len(dlugosci)

if n == 0:
    mediana = 0
elif n % 2 == 1:
    mediana = dlugosci[n // 2]
else:
    mediana = (dlugosci[n // 2 - 1] + dlugosci[n // 2]) / 2

print(liczba_cds)
print(liczba_funkcja_znana)
print(liczba_funkcja_nieznana)
print(mediana)
