from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd
import tkinter as tk
from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk
import os
from tkinter.filedialog import askopenfilename

initialdir=os.getcwd()

def data1():
    global filename1
    filename1 = askopenfilename(initialdir=initialdir,title = "Select target database")
    firstEntryTabThree.delete(0, END)
    firstEntryTabThree.insert(0, filename1)
    global c_
    c_,d_=os.path.splitext(filename1)
    global file1
    file1 = pd.read_csv(filename1)
    global col1
    col1 = list(file1.head(0))
    
def data2():
    global filename2
    filename2 = askopenfilename(initialdir=initialdir,title = "Select query compounds")
    secondEntryTabThree.delete(0, END)
    secondEntryTabThree.insert(0, filename2)
    global file2
    file2 = pd.read_csv(filename2)
    
def file_process(smiles,cpds):
    c_smiles,c_cpd, inv = [],[],[]
    for ds,cpd in zip(smiles,cpds):
        try:
           cs = Chem.CanonSmiles(ds)
           c_smiles.append(cs)
           c_cpd.append(cpd)
        except:
           inv.append(ds)
    return c_smiles, c_cpd, inv

def data_prep(fps1,fps2, t_cpd, q_cpd, t_smiles, q_smiles):
    qu, ta, sim, quc, tac = [], [], [], [], []
    for n in range(len(fps1)): # -1 so the last fp will not be used
       s = DataStructs.BulkTanimotoSimilarity(fps1[n], fps2[:])
       for m in range(len(s)):
           tac.append(t_cpd[n])
           quc.append(q_cpd[:][m])
           ta.append(t_smiles[n])
           qu.append(q_smiles[:][m])
           sim.append(s[m])
    return qu, ta, sim, quc, tac

def submit():
    df1Smiles=file1[file1.iloc[:,0:1].columns[0]]
    df1Cpds=file1[file1.iloc[:,1:2].columns[0]]
    df2Smiles=file2[file2.iloc[:,0:1].columns[0]]
    df2Cpds=file2[file1.iloc[:,1:2].columns[0]]
    t_smiles, t_cpd, t_inv=file_process(df1Smiles,df1Cpds)
    q_smiles, q_cpd, q_inv=file_process(df2Smiles,df2Cpds)
    
    mst = [Chem.MolFromSmiles(x) for x in t_smiles]  ##Convert smiles to drugbank compounds
    msq = [Chem.MolFromSmiles(x) for x in q_smiles] ##Convert smiles to query ligands
    if Criterion2.get()==1:
       fpst = [FingerprintMols.FingerprintMol(x) for x in mst]  #calculare fingerprints of drugbank compounds
       fpsq = [FingerprintMols.FingerprintMol(x) for x in msq]   #calculare fingerprints of query lignads
       chk='rdk'
    elif Criterion2.get()==2:
        fpst = [AllChem.GetMorganFingerprint(x,2) for x in mst]
        fpsq = [AllChem.GetMorganFingerprint(x,2) for x in msq]
        chk='ecpf4'
    elif Criterion2.get()==3:
        fpst = [AllChem.GetMorganFingerprint(x,4) for x in mst]
        fpsq = [AllChem.GetMorganFingerprint(x,4) for x in msq]
        chk='ecpf8'
    elif Criterion2.get()==4:
        fpst = [AllChem.GetMorganFingerprint(x,2,useFeatures=True) for x in mst]
        fpsq = [AllChem.GetMorganFingerprint(x,2,useFeatures=True) for x in msq]
        chk='fcpf4'
    elif Criterion2.get()==5:
        fpst = [AllChem.GetMorganFingerprint(x,4,useFeatures=True) for x in mst]
        fpsq = [AllChem.GetMorganFingerprint(x,4,useFeatures=True) for x in msq]
        chk='fcpf8'
    qu, ta, sim, quc, tac=data_prep(fpst,fpsq, t_cpd, q_cpd, t_smiles, q_smiles)
    d = {'query_cpd': quc, 'query':qu, file1.iloc[:,1:2].columns[0]: tac,'target':ta, 'Similarity':sim}
    df_final = pd.DataFrame(data=d)
    df_final = df_final.sort_values('Similarity', ascending=False)
    
    if Criterion.get()==1:
       df1Act=file1.iloc[:,1:3]
       df_final1=pd.merge(df_final, df1Act, on=file1.iloc[:,1:2].columns[0], how='left')
    elif Criterion.get()==2:
        df_final1=df_final
    df_final1.to_csv(str(c_)+'_'+chk+'_hits.csv', index=False)

form = tk.Tk()
form.title("Similarity Search Tool")
form.geometry("630x230")
tab_parent = ttk.Notebook(form)
tab1 = tk.Frame(tab_parent) #background='#ffffff')
tab_parent.add(tab1, text="SST")


firstLabelTabThree = tk.Label(tab1, text="Select target database file (*.csv)",font=("Helvetica", 12))
firstLabelTabThree.place(x=45,y=10)
firstEntryTabThree = tk.Entry(tab1, width=40)
firstEntryTabThree.place(x=300,y=13)
b3=tk.Button(tab1,text='Browse', command=data1,font=("Helvetica", 10))
b3.place(x=550,y=10) 

secondLabelTabThree = tk.Label(tab1, text="Select query compounds file (*.csv)",font=("Helvetica", 12))
secondLabelTabThree.place(x=45,y=40)
secondEntryTabThree = tk.Entry(tab1,width=40)
secondEntryTabThree.place(x=300,y=43)
b4=tk.Button(tab1,text='Browse', command=data2,font=("Helvetica", 10))
b4.place(x=550,y=40)

Criterion_Label = ttk.Label(tab1, text="Do you have activity column in target molecules?",font=("Helvetica", 12))
Criterion = IntVar()
Criterion.set(2)
Criterion_yes = ttk.Radiobutton(tab1, text='Yes', variable=Criterion, value=1)
Criterion_no = ttk.Radiobutton(tab1, text='No', variable=Criterion, value=2)
Criterion_Label.place(x=30,y=80)
Criterion_yes.place(x=400,y=80)
Criterion_no.place(x=450,y=80)

Criterion_Label2 = ttk.Label(tab1, text="Fingerprint type",font=("Helvetica", 12))
Criterion2 = IntVar()
Criterion2.set(1)
Criterion_rdk = ttk.Radiobutton(tab1, text='RDK', variable=Criterion2, value=1)
Criterion_ecfp4 = ttk.Radiobutton(tab1, text='ECFP4', variable=Criterion2, value=2)
Criterion_ecfp8 = ttk.Radiobutton(tab1, text='ECFP8', variable=Criterion2, value=3)
Criterion_fcfp4 = ttk.Radiobutton(tab1, text='FCFP4', variable=Criterion2, value=4)
Criterion_fcfp8 = ttk.Radiobutton(tab1, text='FCFP8', variable=Criterion2, value=5)
Criterion_Label2.place(x=30,y=120)
Criterion_rdk.place(x=160,y=120)
Criterion_ecfp4.place(x=250,y=120)
Criterion_ecfp8.place(x=350,y=120)
Criterion_fcfp4.place(x=450,y=120)
Criterion_fcfp8.place(x=550,y=120)

b7=Button(tab1, text='Submit', command=submit,bg="orange",font=("Helvetica", 10),anchor=W, justify=LEFT)
b7.place(x=340,y=150)

tab_parent.pack(expand=1, fill='both')


form.mainloop()
