import asyncio
from pipeline import endtoend
from pipeline.meta import BamMeta, ExperimentMeta

ROOT = '/data/encode/H3K27me3/ENCSR793NQA'
TARGERT = "H3K27me3"
ACCESSION = "ENCSR793NQA"

bam = [
    # biological replicates
    BamMeta("breplica1", TARGERT, "ENCFF933CPQ", paired=True, readlen=101),
    BamMeta("breplica2", TARGERT, "ENCFF825SKJ", paired=True, readlen=101),

    # control
    BamMeta("control1", "control", "ENCFF941KVW", paired=True, readlen=101),
    BamMeta("control2", "control", "ENCFF379ZDL", paired=True, readlen=101)
]
bam = {meta.accession: meta for meta in bam}

experiments = [
    ExperimentMeta("breplica1", TARGERT, ACCESSION, [bam["ENCFF933CPQ"]], [bam["ENCFF941KVW"], bam["ENCFF379ZDL"]]),
    ExperimentMeta("breplica2", TARGERT, ACCESSION, [bam["ENCFF825SKJ"]], [bam["ENCFF941KVW"], bam["ENCFF379ZDL"]]),
]
asyncio.run(endtoend.run(ROOT, experiments))

# 1. Делать markdup через sambamba
# 2. paired-end сортировать по имени перед и пулить риды которые парные bedpe
# 3. https://www.biostars.org/p/149119/
# Голь на выдумки хитра. Нужно делать pipelin-ы по типу bash, иначе все совсем грустно становится -
# сожрется вся память, я и не замечу....
# Как тут аккуратно обработать Эти штуковины?...

# ТЗ
# Запуск приложений, в pipelin-e? Если их при этом await, то они полностью выполнятся.
# Как? аннотировать обычные команды о том, что есть файлы-pipes, а что нет.
# В норме процессы просто сразу запускаются, тут же, по-хорошему нужно что-то вроде
# построения, разруливания, запуска. И просто await, если все нормально.
# по хорошему нужно понимать, что есть вход, а что выход. Чтобы просто и понятно разрулить в дальнейшем
# это все несколько печально, на самом деле. Но нет других вариантов, BASH - это звиздец.
# А как их chain при этом? Я же хочу еще и не просто тупые и линейные штуки, мне нужны целые деревья, на самом деле.
# @pipe(a='-', b='-', writes=saveto)
# pipe - штука, которая отвечает за построение дерева? Вроде бы да.
# Вот позвал нашу функцию кто-то, что мы возвращаем? Await[string] - результат.
# И для пользователя теперь можно await и будут чудеса. А под капотом что происходит?
# По капотом... В следующий раз, когда вызывается pipe он видит, что один из параметров - это Await[string], т.е.
# какой-то хитрый объект. pipe нужно сказать, чтобы этот объект обязательно писал в named pip и получить от него конец
# трубы для чтения. Другой pipe в свою очередь смотри куда он там пишет и если это файл, то просто возращает его путь - тогда другой pipe должен догадаться о том, чтобы полностью await
# предыдущий процесс.

# if isawaitable(argument) and argument.ispipeable():
#     argument = argument.__makepipe()  # pipe
# return func(*args, **kwargs)



