flag = input("Enter the flag:")

te = bin(int(flag)) # doit etre compris entre 0 et 4095 => mettre une securite
te = te[2:]
print(te)
print(type(te))
te = list(te)
print(type(te))
print(len(te))
print(te)

if len(te) < 11:
    print("Inferieur")
    ajout = 11 - len(te) # chancer 11 pour 12
    print(ajout)
    print(type(ajout))
    for t in range(ajout):
        te.insert(0,'0')
    print(te)

print("\n")

if 1 == int(te[-1]):
    print("Read paired")
    
if 1 == int(te[-2]):
    print("Read mapped in proper pair")
    
if 1 == int(te[-3]):
    print("Read unmapped")
    
if 1 == int(te[-4]):
    print("Mate unmapped")
    
if 1 == int(te[-5]):
    print("Read reverse strand")
    
if 1 == int(te[-6]):
    print("Mate reverse strand")
    
if 1 == int(te[-7]):
    print("First in pair")
    
if 1 == int(te[-8]):
    print("Second in pair")

if 1 == int(te[-9]):
    print("Not primary alignment")
    
if 1 == int(te[-10]):
    print("Read fails platform/vendor quality checks")
    
if 1 == int(te[-11]):
    print("Read is PCR or optical duplicate")

    # ajouter if 1 == int(te[12]):
