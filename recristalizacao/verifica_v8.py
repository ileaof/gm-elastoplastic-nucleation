# Verificação leve do v8_final.docx: seções 5.7-5.9, figuras e tabelas na ordem certa.
from docx import Document

d = Document(r'C:\Users\ileao\OneDrive\Documentos\Paper Nucleation\Manuscript_ElastoPlastic_Surface_NonEquilibrium_Nucleation_v8_final.docx')
achados = {}
ordem = []
for i, p in enumerate(d.paragraphs):
    txt = p.text.strip()
    if not txt:
        continue
    for alvo in ('5.7.', '5.8.', '5.9.', 'Figure 7.', 'Figure 8.', 'Figure 9.', 'Figure 10.',
                 'Figure 11.', 'Table 8.', 'Table 9.', 'Table 10.', '6. Conclusions'):
        if txt.startswith(alvo):
            achados[alvo] = i
            ordem.append((i, alvo, txt[:95]))
print(f"parágrafos: {len(d.paragraphs)}; imagens inline: {len(d.inline_shapes)}")
for i, alvo, txt in ordem:
    print(f"{i:5d}  {alvo:11s} {txt}")
assert achados['5.9.'] > achados['5.8.'] and achados['Figure 11.'] > achados['Figure 10.'], 'ordem errada'
print('OK: 5.9 e Figure 11 depois de 5.8/Figure 10 e antes de', achados['6. Conclusions'])