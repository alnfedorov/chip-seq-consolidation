## What is the pain the pain to address?
How to reliably merge several biological replicates? Especially if they don't have a lot of reads.
* IDR - only two replicates in practice, bad for broad histon marks
* Pool - doesn't account for biases in the specific platforms experiments. Might be good in some cases.
* Peak-calling and merge later. Good option, but how reliable is it? How to merge?

For ULI CHip-seq it is common to have few reads for the peak-source histone marks. Hard to reliably peak-call ~10m reads. 

CHip-seq has it's own problems - antibody quality, sequencing and alignment problems. Errors on each step introduces 
unnecessary biases that affect peak-calling algorithms.

What do I want? An option to reliably merge several biological replicates. Even if there is problems during peak-calling 
or any other pipeline steps. How to do this? Peaks MUST have some sort of biological similarity inside/around them in order to 
distinguishable from other genomic locations. Hence we can relay on the ability of neural networks to learn in presence 
of the noisy labels to filter out false-positive and probably false-negative regions.
## QC options
### Tuning sandbox
##### Simulations:  
Simulate chip-seq experiment from scratch.  
* Requires high-quality gt regions(not available for histone's chip-seq).  
* Sophisticated and perhaps unreliable tooling.  
##### Real data:
Use ENCODE data.  
* Low-quality and noisy.  
* There is no gt regions. Hence, it is hard to measure the quality of the reconstruction.  
##### Hybrid approach
Use deeply sequenced epigenomes as a golden standart and apply several peak-calling tools to obtain reliable set of 
peaks. Use them as a ground true for simulations.  
* More or less reliable approach.  
* Use simulations with random set of peaks as a check for adequacy in simulations.  
### Metrics
1. IOU between replicates and consensus.
2. Precision and recall for sandbox replicas.
3. Fraction of reads in peaks?
4. Signal to noise ratio?
5. Correlation with RNA-seq and active genes
# TODO
1. Надо прочитать chipulate, чтобы лучше понимать проблемы симуляций.
2. Надо запилить еще и CHIP-S или как-то так и посмотреть на его результат.
3. Надо понять ЧТО ДАНО, ЧТО ЕСТЬ ВЫХОД. Это самое важное.


#### ИТОГО, балансировка участков по:
1. Доля обогащения
2. СКО отношения числа ридов (меня не слишком интересуют константные участки)
3. Геномная разметка (экзон/интрон/UTR и т.д.)
4. Культура
5. Метка
6. Число уникальных ридов

Для каждой комбинации классов заводим лист, туда складываем индексы соответствующих элементов. 
Потом равновероятно выбираем комбинацию и рандомно сэмплируем индекс объекта от туда.

Как хранить сами объекты? Очень хороший вопрос.  
Список регионов с геномной разметко и пересечением с консенсусом.  
Для каждой метки есть список экспериментов где у каждого эксперимента есть список регионов идентичный первому.


Список регионов -> (мета региона, метка -> (мета эксперимента, отсортированные регионы))



#### АУГМЕНТАЦИЯ.
1. Движение окна на +- length / 2
2. Симулирование доп. ридов из контроля или экспериментов. ИМЕННО СИМУЛИРОВАНИЕ, не тупой +рандом
3. Вроде бы все.

Разбивать геном будем на не пересекающиеся бины по размеру окна. Так значительно проще.   
И рандомно сэмплируем эти окна. Собственно, разделение ответственности понятно. Датасет будем делать через 
каналы, следовательно, все будет достаточно просто.  


Для каждого окна считаем статистику. Затем все значения статистик нормализуем, разбивая на квартили   
Декартово произведение всех классов. Если какой-то из классов не встречается, то перенормализуем.  


Ок, давай сделаем красиво!


What to i need to do?
1. Нужно сделать функцию, которая все это процессит в такие вот красивые файлики.
2. Нужен пайплайн, который все это сделает. Т.е. на входе ID из encode,
а на выходе результаты peak-calling и bigWig c ридами в контроле/ридами вне контроля
Какое окно? Пусть 1001bp * 25 = 25kbp
Ок, но с окном есть проблема. Как тогда предсказывать риды из такого окна? Большая проблема
Все фигня. Пока что мы пытаемся выбрать окно. Или не выбрать.

Как построить процесс выборки правильным образом?
По-хорошему мы имеем дело с одним большим изображением.
Какие могут быть соображения?
1. Число примеров обогащение / не обогащение должно быть сопоставимо в батче.
2. Примеры с обогащением должны быть сложными.
   Т.е. должны быть в том числе ситуации, когда обогащения мало - обошащения много и т.д.
3. Константные участки меня не интересуют вообще. От слова совсем. Максимум парочку добавить, не больше.
4. Мне дано окно - размер входа модели. Размер окна фиксирован.
5. В идеале мне нужно для каждой позиции в геноме посчитать есть ли там обогащенные участки и от этого плясать.
   Получается, что для каждой позиции я имею float - вероятность выбора позиции в тех или иных обстоятельствах.
6. По мимо этого меня интересуют разное ИСХОДНОЕ количество ридов.


Отношение ридов - это классно, конечно. И нам интересно именно оно. Но вот тулам интересны абсолютные величины.
Что должна делать модель? Она должна аккуратно регрессировать правильное отношение между ридами контроля и эксп-а в данной области.
т.е. A/B -> A'/B'. Что тогда получаем? Мы знаем скорректированное отношение. И пусть мы знаем сумму всех A'_i и B'_i
min (A''_i - B''_i) - (A'_i / B'_i) given sum{A''_i} is bound and sum{B''_i} is bound and A''_i > 0 and B''_i > 0.
На самом деле это не сложная задача, достаточно взять множители лагранжа и посмотреть, что получится.
Нет, не достаточно. Линейные оптимизации могут дать любой хлам т.к. мы аппроксимируем отношение.
Нужно рассматривать что-то вроде пошаговой оптимизации базового отношения X_i/Y_i и новых ридов на каком-то шаге
Т.е. X_i + a / Y_i + b

Можно ли как-то аппроксимировать эту штуку? Хороший вопрос. Очень.
Но нас для работы интересуют именно абсолютные величины - MACS2 peak calling работает именно на них.

signal track construction
1. macs2 callpeak -t H3K36me1_EE_rep1.bam -c Input_EE_rep1.bam  -B --SPMR -g hs -n H3K36me1_EE_rep1
2. macs2 bdgcmp -t H3K36me1_EE_rep1_treat_pileup.bdg -c H3K36me1_EE_rep1_control_lambda.bdg -o H3K36me1_EE_rep1_FE.bdg -m FE

Подход не очень, лучше делать через genomecov? И быстрее и лучше
BAM должен быть отсортирован по позиции !!!!
genomecov -bga -ibam file.bam > breplica1.bdg
bedtools slop -i breplica1.bdg -g chromInfo.txt -b 0 > breplica1.slop.bdg
bedClip breplica1.slop.bdg chromInfo.txt breplica1.clip.slop.bdg
LC_COLLATE=C sort -k1,1 -k2,2n ${F}.clip > ${F}.sort.clip             # sort by something
bedGraphToBigWig ${F}.sort.clip ${G} ${F/bdg/bw}