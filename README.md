В данной работе использованы данные: транскриптом дрожжей Saccharomyces cerevisiae - реакция на 45 мМ L-лактата (https://doi.org/10.1007/s00253-023-12863-z)

# Работа с удаленным сервером bash

Для работы использована программа PuTTY

- ssh -p 627 2025_RR_St_X@bioinformatics.isu.ru

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


Скачивание данных для дрожжей: 

```bash
fasterq-dump --threads 2 -A --progress SRR24466389; fasterq-dump --threads 2 -A --progress SRR24466390; fasterq-dump --threads 2
-A --progress SRR24466391; fasterq-dump --threads 2 -A --progress SRR24466374; fasterq-dump --threads 2 -A --progress
SRR24466375; fasterq-dump --threads 2 -A --progress SRR24466376
```

## Выравнивание чтений на референс: дрожжи

### download reference / скачиваем референс

```bash
wget https://ftp.ensembl.org/pub/release-108/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz
Wget
https://ftp.ensembl.org/pub/release-108/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.tople
vel.fa.gz
```

### unpack archive / распаковываем архивы

```bash
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz
```

### hisat: build index and prepare splice file / создаём индекс и готовим файл с данными сплайсинга в hisat2

```bash
hisat2-build Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa yeast_index
hisat2_extract_splice_sites.py Saccharomyces_cerevisiae.R64-1-1.108.gtf > yeast_splice_sites.txt
```

### hisat: align & samtools: make sorted bam / выравниваем с помощью hisat2 и сортируем bam-файл с помощью samtools

```bash
for sample in `ls *_1.fastq`; do base=$(basename $sample "_1.fastq"); hisat2 -x yeast_index --known-splicesite-infile yeast_splice_sites.txt -p 8 -1 ${base}_1.fastq -2 ${base}_2.fastq | samtools view --threads 2 -bS | samtools sort --threads 2 -o $base.bam; done
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

