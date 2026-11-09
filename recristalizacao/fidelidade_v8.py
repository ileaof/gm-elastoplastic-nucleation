# Fidelidade do v8 vs v7: a únicas diferenças de texto devem ser o bloco 5.9 (Fig 11 + Tabela 10);
# e os números da Tabela 10 devem ser os da análise.
import sys
sys.stdout.reconfigure(encoding='utf-8')
from docx import Document

V7 = r'C:\Users\ileao\OneDrive\Documentos\Paper Nucleation\Manuscript_ElastoPlastic_Surface_NonEquilibrium_Nucleation_v8_final.docx'
V8 = r'C:\Users\ileao\OneDrive\Documentos\Paper Nucleation\Manuscript_ElastoPlastic_Surface_NonEquilibrium_Nucleation_v7_total.docx'
t7 = [p.text for p in Document(V7).paragraphs]
t8 = [p.text for p in Document(V8).paragraphs]
print(f"v7: {len(t7)} parágrafos; v8: {len(t8)}")
d8 = [x for x in t8]
for p7 in t7:
    if p7 in d8:
        d8.remove(p7)
print("--- em v8 e não em v7 (esperado: só o bloco 5.9) ---")
for x in d8:
    print(" •", x[:100].replace('\n', ' '))
for x in t7:
    if None or (x and x not in t8):
        print(" REMOVIDO do v7:", x[:100])

doc8 = Document(V8)
print("\n--- Tabela 10 (última tabela do v8) ---")
tb = doc8.tables[-1]
for linha in tb.rows:
    print(" | ".join(c.text for c in linha.cells))
print("\n--- Table 9: última linha + contagem ---")
for t9 in doc8.tables[-2:]:
    pass
tb9 = doc8.tables[-2]
print("linhas:", len(tb9.rows), "; última:", " | ".join(c.text for c in tb9.rows[-1].cells))
print("\n--- Table 8: restos armazenados (coluna 5) ---")
tb8 = doc8.tables[-3]
for linha in tb8.rows[1:]:
    print(" | ".join(c.text for c in linha.cells))