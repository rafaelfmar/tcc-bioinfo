## Gerar e atualizar dados GO e KEGG
Antes de processar os dados é necessário realizar o download dos arquivos necessários.<br/>
No windows:
```
UpdateGo.ps1
UpdateKo.ps1
```
No Linux:
```
UpdateGo.sh
UpdateKo.sh
```

Após o download deve-se executar os scripts:<br/>
Processar dados GO:
```
dag.R
```
Processar dados KEGG:
```
ortholog.R
```

## Calcular a similaridade semântica entre dois termos
Para calcular a similaridade entre termos deve-se executar o script informando o nome de ambos os termos.<br/>
Comando:
```
Rscript ./termComb.R <termo A> <termo B>
```
Exemplo:
```
Rscript ./termComb.R GO:0005975 GO:1901135
```
## Calcular a similaridade semântica de um ou mais pares de genes
Para calcular a similaridade entre um ou mais pares de genes deve-se executar o script informando o caminho para o arquivo CSV contendo os termos.<br/>
Comando:<br/>
```
Rscript ./genePairComb.R <arquivo csv contendo lista de genes com seus respectivos termos>
```
Exemplo:<br/>
```
Rscript ./genePairComb.R ./data/input/E.coli/assimilatory_sulfate_reduction.csv
```
O resultado será exibido no terminal, e também armazenado na pasta ./data/result/.<br/>
## Calcular a similaridade entre genes e classificá-los de acordo com o resultado
Para calcular a similaridade entre genes e e classificá-los deve-se executar o script informando o caminho para a pasta contendo os arquivos CSV com os termos e a pasta em que os resultados serão armazenados.<br/>
Comando<br/>
```
Rscript ./geneListComb.R <pasta contendo arquivos csv> <pasta que será armazenado o resultado de classificação>
```
Exemplo:<br/>
```
Rscript ./geneListComb.R ./data/input/E.coli ./data/result/E.coli
```
