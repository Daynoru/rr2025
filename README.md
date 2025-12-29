# Резюме

*Saccharomyces cerevisiae* (пекарские дрожжи) — модельный эукариотический организм, используемый для фундаментальных исследований и производства широкого спектра соединений. Несмотря на обширные данные о реакции дрожжей на различные стрессоры, систематический анализ их транскрипционного ответа на L-лактат — соединение, представляющее как промышленный, так и физиологический интерес, — до сих пор отсутствовал.

Данный репозиторий содержит материалы исследования, в котором впервые проведён детальный анализ изменений в экспрессии генов *S. cerevisiae* под воздействием 45 мМ L-лактата. Основной целью работы было восполнение пробела в данных и выявление специфических черт клеточного ответа на данное соединение.

В репозитории представлена последовательность шагов (pipeline) для биоинформатического анализа данных RNA-seq, нацеленного на выявление дифференциальной экспрессии генов. Он включает в себя руководство по использованию командной оболочки Bash, выравниванию (HISAT2), количественной оценке экспрессии (featureCounts) и статистическому анализу (пакет DESeq2 в R). Этот конвейер обеспечивает воспроизводимую основу для обработки транскриптомных данных, аналогичных полученным в данном исследовании.

<p align="center">
  <img src="https://github.com/Daynoru/rr2025/blob/main/image/Full_workflow.png" alt="workflow image" width="70%"/>
</p>


# Материалы и методы

Для работы с удаленным сервером bash была использована программа [PuTTY версии 0.83](https://www.chiark.greenend.org.uk/~sgtatham/putty/releases/0.83.html). Основной анализ данных выполнен на языке программирования [R версии 4.3.1](https://www.r-project.org/)  в программной среде [RStudio](https://posit.co/download/rstudio-desktop/). Полученная таблица подсчётов была обработана с использованием пакета [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) версии 1.50.2 с целью сравнения уровней экспрессии. Визуализация данных выполнена с помощью пакетов [ggplot2](https://ggplot2.tidyverse.org/) версии 4.0.1, [enhancedVolcano](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html) версии 1.28.2 и [DEGReport](https://www.bioconductor.org/packages/release/bioc/html/DEGreport.html)  версии 1.46.0 для R. Для скачивания данных с сервера использована программа [Filezilla](https://filezilla-project.org/) версии 1.12.1.

Референсные последовательности доступны по [данной ссылке](https://ftp.ensembl.org/pub/release-108/gtf/saccharomyces_cerevisiae/). Данные секвенирования доступны под номером доступа [PRJNA970274](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA970274&o=acc_s%3Aa). Использованные в работе данные ответа на 45 мМ L-лактата доступны под номерами: SRR24466389; SRR24466390; SRR24466391; SRR24466374; SRR24466375; SRR24466376.

# Скачивание и выравнивание чтений на референс

Перед непосредственным анализом дифференциальной экспрессии генов было необходимо произвести скачиваение и выравнивание чтений на референсную последовательность. Для данной работы был использован удаленный сервер bash. 

Bash (Bourne-Again SHell) — это командная оболочка и язык сценариев для Linux, macOS и других UNIX-подобных систем. Позволяет пользователю взаимодействовать с системой через команды, автоматизировать задачи и писать скрипты для управления процессами, файлами и программами.

## Работа с удаленным сервером bash

Для входа в аккаунт на сервере используем команду:

```bash
- ssh -p номер_порта имя_пользователя@адрес_сервера
```

### Основные команды bash

| Команда | Смысл | Опции и комментарии |
|-------------|-------------|-------------|
| pwd   | Где я? (print working directory)    |  |
| cd | Перейти в другое расположение (change directory)    |. = начало относительного пути; .. = подняться на уровень выше; ~ = домашняя папка; / = разделитель уровней; Tab = автозаполнение      |
| mkdir   | Создать пустую папку (make directory)   | |
| ls     | Посмотреть содержимое папки (ls -lh /home/user/test)  | -h удобное отображение размера; -t сортировка по размеру; -l табличкой  |
| head, tail     | Посмотреть начало или конец файла   |  |
| cat | Посмотреть весь файл   |  |
| grep  | Поиск текста в файле|   |
| >  | Записать вывод команды в файл |  > file.txt |
| mv, cp  | Переместить или переименовать; копировать | -r = работа с папкой  |
| rm  | Удалить| -r = удаляет папку  |
| nano  | Текстовый редактор|  |
| sudo название_команды опции файл  | Права администратора |  |
| scp /path/to/local/file username@server:/path/to/file | Скачивание между локальным компьютером и удалённым сервером |  |
| unzip file.zip | Архивы |  |
| screen -S myprocess | Создание screen'а | Ctrl+A+D #чтобы выйти; screen -r myprocess #вернуться; screen -ls #показать список |

Скачивание данных для дальнейшей работы было осуществлено с использованием инструмента **fasterq-dump**, который позволяет скачивать данные из базы данных NCBI's Sequence Read Archive (SRA): 

```bash
fasterq-dump --threads 2 -A --progress SRR24466389; fasterq-dump --threads 2 -A --progress SRR24466390; fasterq-dump --threads 2
-A --progress SRR24466391; fasterq-dump --threads 2 -A --progress SRR24466374; fasterq-dump --threads 2 -A --progress
SRR24466375; fasterq-dump --threads 2 -A --progress SRR24466376
```

Скачивание референсных последовательностей:

```bash
wget https://ftp.ensembl.org/pub/release-108/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz
Wget
https://ftp.ensembl.org/pub/release-108/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.tople
vel.fa.gz
```

## Выравнивание чтений на референс

Архивы с референсным геномом и аннотациями были распакованы с помощью команды **gunzip**. Полученные файлы были использованы для построения индекса и подготовки данных о сайтах сплайсинга для программы **hisat2**.

```bash
#распаковка референсного генома и аннотаций Saccharomyces cerevisiae
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz

#построение индекса референсного генома для hisat2
hisat2-build Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa yeast_index

#извлечение информации о сайтах сплайсинга из аннотаций
hisat2_extract_splice_sites.py Saccharomyces_cerevisiae.R64-1-1.108.gtf > yeast_splice_sites.txt
```

Выравнивание RNA-seq ридов на референсный геном было выполнено с помощью **hisat2** с последующей конвертацией и сортировкой результатов в формат BAM посредством **samtools**.

```bash
for sample in `ls *_1.fastq`; do 
    base=$(basename $sample "_1.fastq")
    
    # -x: индекс референсного генома
    # --known-splicesite-infile: файл с известными сайтами сплайсинга
    # -p 8: использование 8 потоков для ускорения вычислений
    # -1 и -2: файлы с прямыми и обратными ридами
    hisat2 -x yeast_index --known-splicesite-infile yeast_splice_sites.txt -p 8 \
           -1 ${base}_1.fastq -2 ${base}_2.fastq | \
    
    #конвертация выравниваний из формата SAM в BAM
    samtools view --threads 2 -bS | \
    
    #сортировка выравниваний по координатам в геноме
    samtools sort --threads 2 -o ${base}.bam
done
```

# Анализ дифференциальной экспрессии генов

 Для количественной оценки экспрессии генов была использована программа **featureCounts**. Подсчёт ридов, выравненных на аннотированные гены, выполнялся с указанием определенных параметров:

```bash

export PATH=$PATH:/media/secondary/apps/subread-2.0.4-Linux-x86_64/bin


# -s 2: риды считываются в обратной ориентации (reverse strand)
# -T 2: использование 2 потоков (threads) для ускорения вычислений
# -p: учёт парных ридов (paired-end reads)
# -a: файл аннотаций в формате GTF
# -o: имя выходного файла с результатами
# $(ls *.bam): все BAM-файлы в текущей директории как входные данные
featureCounts -s 2 -T 2 -p -a Saccharomyces_cerevisiae.R64-1-1.108.gtf \
-o allSamples.featureCounts.txt $(ls *.bam)
```
Для скачивания полученных данных использовали программное обеспечение **Filezilla**.


## Основные команды базового R
```R
setwd() #установка рабочей директории
install.packages() #установка пакетов из CRAN
library() #загрузка установленного пакета в сессию
str() #структура объекта
unique() #уникальные значения
hist() #построение гистограммы
mean() #среднее арифметическое
median() #медиана
shapiro.test() #тест на нормальность
wilcox.test() #непараметрический тест для сравнений
```
### Команды пакетов для R
Используемые библиотеки: openxlsx, ggplot2, ggpubr, scales.
```R
read.xlsx() #чтение данных из Excel файлов
ggplot() #создание графиков грамматикой графиков
ggsave() #сохранение ggplot графиков в файлы
```

На  основе данных о влиянии кислотной (D-молочная кислота) и нейтральной (D-лактат натрия) форм D-изомера на рост дрожжей и экспрессию ключевых генов был построен график:
```R
yeast.growth <- read.xlsx("253_2023_12863_MOESM2_ESM.xlsx", sheet = 3)
ggplot(yeast.growth, aes(x=Condition, y=RGR2)) + #данные для графика
  geom_boxplot(color = "forestgreen", fill = "pink") + #отрисовка боксплота
  geom_pwc(method="wilcox_test", label = " {p.adj}", 
            p.adjust.method = "holm") + #статистические данные
geom_jitter(width = 0.1, colour = "forestgreen") 
    

ggsave("results1.png", device=png, width=16, height=12, units="cm") #сохранение результатов 
```

<p align="center">
<img src="https://github.com/Daynoru/rr2025/blob/main/image/results1.png" alt="График" width="50%" />
</p>

### Визуализация дифференциальной экспрессии генов

Идентификация дифференциально экспрессируемых генов (DEGs) между образцами, обработанными L-лактатом, и контрольными образцами были выполнены в среде R с использованием пакета DESeq2. Визуализация результатов осуществлена с помощью пакета EnhancedVolcano.

```R
library(BiocManager)

setwd("D:\\uni\\RR6") #установка рабочей директории
count_table <- read.delim("allSamples.featureCounts.txt", skip=1, row.names="Geneid") #загрузка и подготовка данных для анализа
sample_table <- data.frame(condition=c("DL", "DL", "DL", "control",
                                       "control"))
library(DESeq2)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count_table[,6:10], colData = sample_table, design = ~ condition) #объединяем данные подсчётов с метаданными образцов
dds <- DESeq(ddsFullCountTable) #основная функция анализа: нормализация и статистический тест
res <- results(dds)

library(EnhancedVolcano)
EnhancedVolcano(res, lab = rownames(res),
                x = 'log2FoldChange', y = 'pvalue',
                pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                title="Large Title", subtitle="Subtitle",
                col = c("grey30", "grey30", "grey30", "red2"),
                xlab="", ylab = bquote(~-Log[10] ~ italic(p)),
                caption="", selectLab = "", legendPosition = 'none')

DEGs <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05 & complete.cases(res$padj), ] #фильтрация и сохранение значимо изменённых генов
DEGs <- DEGs[order(DEGs$log2FoldChange), ] #сортировка по величине изменения

library(openxlsx)
DEGs$Transcript <- row.names(DEGs)
write.xlsx(x = DEGs, file = "DEGs_yeast.xlsx") #сохранение в файл Excel

EnhancedVolcano(res, lab = rownames(res), #визуализация результатов - вулканоплот
                x = 'log2FoldChange', y = 'pvalue',
                pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                title="S. cerevisiae", subtitle="45 mM L-lactate vs control",
                col = c("grey30", "grey30", "grey30", "red2"),
                xlab="", ylab = bquote(~-Log[10] ~ italic(p)),
                caption="", 
                selectLab = rownames(DEGs),  #выделение на графике всех значимых генов
                legendPosition = 'none')
```


В результате анализа были идентифицированы дифференциально экспрессируемые гены. На вулканоплоте показано распределение всех генов по величине изменения экспрессии (log2FoldChange, ось X) и статистической значимости (-log10(p-value), ось Y). Гены, удовлетворяющие критериям значимости (|log2FC| > 1, padj < 0.05), выделены красным цветом.

<p align="center">
<img src="https://github.com/Daynoru/rr2025/blob/main/image/diff_expression_withlabels.png" alt="График" width="50%" />
</p>


Для интерпретации биологического значения выявленных изменений был проведен функциональный анализ с использованием базы данных  [Saccharomyces genome database](https://yeastgenome.org/). Наиболее интересные гены, демонстрирующие значительные изменения экспрессии в ответ на обработку L-лактатом, представлены в таблице.

| Ген | Описание | Пояснение | Изменение |
|-------------|-------------|-------------|-------------|
| [ECM3 / YOR092W](https://yeastgenome.org/locus/S000005618)   | involved in signal transduction and the genotoxic response; induced rapidly in response to treatment with 8-methoxypsoralen and UVA irradiation; relocalizes from ER to cytoplasm upon DNA replication stress    | Снижение экспрессии гена в ответ на повреждение ДНК или окислительный стресс для экономии клеточных ресурсов  |↓ |
| [LEU1 / YGL009C](https://yeastgenome.org/locus/S000002977)   | Isopropylmalate isomerase; catalyzes the second step in the leucine biosynthesis pathway   | Связывает кислотный стресс с регуляцией метаболизма аминокислот. Возможная компенсаторная активация. |↑ |
| [ARN1 / YHL040C](https://yeastgenome.org/locus/S000001032) | ARN family transporter for siderophore-iron chelates; responsible for uptake of iron bound to ferrirubin, ferrirhodin, and related siderophores; protein increases in abundance and relocalizes to the vacuole upon DNA replication stress  | Регуляция трансляции в условиях стресса  |↓ |
|[ZPS1 / YOL154W](https://yeastgenome.org/locus/S000005514)| Putative GPI-anchored protein; transcription is induced under low-zinc conditions, as mediated by the Zap1p transcription factor, and at alkaline pH   |Участвуют в ремоделировании клеточной стенки, усилении барьера или передаче сигналов. Это может быть частью структурной адаптации к стрессу |↑ |
| [snR84](https://yeastgenome.org/locus/S000028466) | H/ACA box small nucleolar RNA (snoRNA); guides pseudouridylation of large subunit (LSU) rRNA at position U2266; overexpression confers resistance to baking-associated stress   | Снижение активности биогенеза/модификации рибосом, что имеет смысл в стрессовых условиях для снижения нагрузки на клетку |↓ |
| [SOD1 / YJR104C](https://yeastgenome.org/locus/S000003865) | Cytosolic copper-zinc superoxide dismutase; also sulfide oxidase; detoxifies superoxide and hydrogen sulfide; stabilizes Yck1p and Yck2p kinases in glucose to repress respiration; phosphorylated by Dun1p, enters nucleus under oxidative stress to promote transcription of stress response genes; abundance increases under DNA replication stress  | SOD1 превращает два супероксид-аниона в менее опасную перекись водорода (H₂O₂) и молекулярный кислород. Это первая и важнейшая линия защиты от окислительного стресса. |↑ |
| [RMD5 / YDR255C](https://yeastgenome.org/locus/S000002663) | Component of GID Complex that confers ubiquitin ligase (U3) activity; necessary for polyubiquitination and degradation of the gluconeogenic enzyme fructose-1,6-bisphosphatase | Вызывает деградацию одного из ключевых ферментов глюконееогенеза, снижая энергозатраты. |↑ |

# Выводы

Таким образом, на основании изменения экспрессии генов *S. cerevisiae* в ответ на 45 мМ L-лактата, мы можем наблюдать, что ответ затрагивает множество метаболических путей. Можно сделать выводы о том, что для дрожжей характерна многоуровневая адаптация к стрессовым условиям, которая направлена как на нейтрализацию непосредственной угрозы, экономию ресурсов, так и перестройку метаболизма всей клетки для выживания. 

Изменение экспрессии проанализированных генов отражает не хаотичную реакцию на стресс, а сложную, многоуровневую и внутренне согласованную адаптивную реакцию. Эта реакция направлена на: 1) нейтрализацию непосредственной угрозы, 2) глобальную экономию ресурсов и 3) перестройку метаболизма для выживания в ущерб росту.
