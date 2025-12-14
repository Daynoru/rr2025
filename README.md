Модельный дрожжевой организм *Saccharomyces cerevisiae* является широко используемым объектом как фундаментальных, так и прикладных исследований, включая разработку биосенсоров и промышленное производство фармацевтических соединений. Однако, несмотря на многочисленные исследования транскрипционного ответа *S. cerevisiae* на различные вещества, реакция на некоторые соединения, производимые в дрожжах, в частности на L-лактат, остаётся неизученной. В данной работе мы исследуем [транскрипционный ответ штамма BY4742](https://doi.org/10.1007/s00253-023-12863-z) на 45 мМ L-лактата.

Данный репозиторий содержит данные и программный код, необходимые для анализа транскрипционного ответа *S. cerevisiae*. 

# Материалы и методы

Для работы с удаленным сервером bash была использована программа [PuTTY release (0.83)](https://www.chiark.greenend.org.uk/~sgtatham/putty/releases/0.83.html). Основной анализ данных выполнен на языке программирования [R v4.3.1](https://www.r-project.org/)  в программной среде [RStudio](https://posit.co/download/rstudio-desktop/). Полученная таблица подсчётов была обработана с использованием пакета DESeq2 (Love и др., 2014) версии 1.34.0 для R (R Core Team, 2022) версии 4.1.2 с целью сравнения уровней экспрессии. Визуализация данных выполнена с помощью пакетов ggplot2 (Wickham, 2016) версии 3.4.2, enhancedVolcano (Blighe и др., 2023) версии 1.12.0 и DEGReport (Pantano и др., 2023) версии 1.30.3 для R.

# Работа с удаленным сервером bash

Для входа в аккаунт на сервере используем команду

```bash
- ssh -p 627 2025_RR_St_X@bioinformatics.isu.ru
```
где X - номер студента. 

### Основные команды

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

Далее представлен код, который был использован в процессе работы с bash.

Скачивание данных для дальнейшей работы: 

```bash
fasterq-dump --threads 2 -A --progress SRR24466389; fasterq-dump --threads 2 -A --progress SRR24466390; fasterq-dump --threads 2
-A --progress SRR24466391; fasterq-dump --threads 2 -A --progress SRR24466374; fasterq-dump --threads 2 -A --progress
SRR24466375; fasterq-dump --threads 2 -A --progress SRR24466376
```

## Выравнивание чтений на референс

### Скачивание референсных последовательностей

```bash
wget https://ftp.ensembl.org/pub/release-108/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz
Wget
https://ftp.ensembl.org/pub/release-108/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.tople
vel.fa.gz
```

### Распаковка архивов

```bash
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz
```

### Создание индекса hisat: build index and prepare splice file / создаём индекс и готовим файл с данными сплайсинга в hisat2

```bash
hisat2-build Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa yeast_index
hisat2_extract_splice_sites.py Saccharomyces_cerevisiae.R64-1-1.108.gtf > yeast_splice_sites.txt
```

### hisat: align & samtools: make sorted bam / выравниваем с помощью hisat2 и сортируем bam-файл с помощью samtools

```bash
for sample in `ls *_1.fastq`; do base=$(basename $sample "_1.fastq"); hisat2 -x yeast_index --known-splicesite-infile yeast_splice_sites.txt -p 8 -1 ${base}_1.fastq -2 ${base}_2.fastq | samtools view --threads 2 -bS | samtools sort --threads 2 -o $base.bam; done
```

### Код для анализа дифференциальной экспрессии генов

```bash
# export
export PATH=$PATH:/media/secondary/apps/subread-2.0.4-Linux-x86_64/bin
## featureCounts: quantify
featureCounts -s 2 -T 2 -p -a Saccharomyces_cerevisiae.R64-1-1.108.gtf \
-o allSamples.featureCounts.txt $(ls *.bam)
```

Пример скачивания данных с сервера. В данной работе для скачивания использовалось программное обеспечение Filezilla

```bash
.\pscp -P 627 2025_RR_StX@bioinformatics.isu.ru:~/allSamples.featureCounts.txt 
```

# Программирование на R

## Программный код

### Основные команды базового R
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
### Команды пакетов
Используемые библиотеки: openxlsx, ggplot2, ggpubr, scales
```R
read.xlsx() #чтение данных из Excel файлов
ggplot() #создание графиков грамматикой графиков
ggsave() #сохранение ggplot графиков в файлы
```

### Пример кода для ggplot2
```R
yeast.growth <- read.xlsx("253_2023_12863_MOESM2_ESM.xlsx", sheet = 3)
ggplot(yeast.growth, aes(x=Condition, y=RGR2)) + #данные для графика
  geom_boxplot(color = "forestgreen", fill = "pink") + #рисуем боксплот
  geom_pwc(method="wilcox_test", label = " {p.adj}", 
            p.adjust.method = "holm") +
geom_jitter(width = 0.1, colour = "forestgreen") 
    

ggsave("results1.png", device=png, width=16, height=12, units="cm")
```
Результат кода:

<img src="https://github.com/Daynoru/rr2025/blob/main/image/results1.png" alt="График" width="50%" />

### Визуализация дифференциальной экспрессии генов

```R
library(BiocManager)

setwd("D:\\uni\\RR6")
count_table <- read.delim("allSamples.featureCounts.txt", skip=1, row.names="Geneid")
sample_table <- data.frame(condition=c("DL", "DL", "DL", "control",
                                       "control"))
library(DESeq2)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count_table[,6:10], colData = sample_table, design = ~ condition)
dds <- DESeq(ddsFullCountTable)
res <- results(dds)

library(EnhancedVolcano)
EnhancedVolcano(res, lab = rownames(res),
                x = 'log2FoldChange', y = 'pvalue',
                pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                title="Large Title", subtitle="Subtitle",
                col = c("grey30", "grey30", "grey30", "red2"),
                xlab="", ylab = bquote(~-Log[10] ~ italic(p)),
                caption="", selectLab = "", legendPosition = 'none')

DEGs <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05 & complete.cases(res$padj), ]
DEGs <- DEGs[order(DEGs$log2FoldChange), ]
library(openxlsx)
DEGs$Transcript <- row.names(DEGs)
write.xlsx(x = DEGs, file = "DEGs_yeast.xlsx")
```

Визуализация результатов анализа:

<img src="https://github.com/Daynoru/rr2025/blob/main/image/diff.png" alt="График" width="50%" />

