import glob

fasta_files = glob.glob("data/*/*.fasta",recursive=True)
fasta_files += glob.glob("data/**/*.fastq",recursive=True)

samples = ["WT","Mut1","Mut2","Mut3","Mut4","Mut5"]
# reads split reads in 9 subsets to manage ~1000x depth sequencing
n_splits = 9
split_numbers = list(range(1,n_splits+1))

# samples=["WT"]
# split_numbers=[9]




# remove duplicate lines that start with the same 4 columns
awk '!seen[$1$2$3$4$5]++' {input} > {output}

# string between 3rd and 4th occurence of "|" in each row in rows not starting with #
 {input} | awk -F "|" '{print $3}' | sort | uniq -c | sort -nr | head -n 10
awk '!/^#/ {print $0}' OkinawaSeq/data/processed/vcf/hifiasm/Mut1.part_001.vs.WT.part_002.ann.vcf | awk -F "|" '{print $4, $10}'

     81 corA c.178dupG
     81 msrA c.603T>C
     63 rcnA c.-61A>T
     81 rhaA c.339G>T
     18 thiM c.-555T>A
     81 corA c.178dupG
     63 rcnA c.-61A>T
     81 rhaA c.339G>T
     18 thiM c.-555T>A
     81 corA c.607_720del
     72 mgtA c.2444C>T
     81 phoP c.256G>C
     63 rcnA c.-10A>G
     18 thiM c.-606T>C
      9 btuB c.901delT
     81 corA c.178dupG
     63 rcnA c.-61A>T
     81 rhaA c.339G>T
     18 thiM c.-555T>A
     81 corA c.185C>A
     81 mgtA c.2336T>G

"""
corA 
185C>A

corA 
178dupG

corA 
607_720del

msrA 
603T>C

rcnA 
-61A>T

  rcnR 
-60T>A

  rcnA 
-10A>G

  rcnR 
-111T>C

  paeA 
-359T>C

  rhaA 
339G>T

  mgtA 
2444C>T

  mgtA 
2336T>G

  phoP 
256G>C

"""

"data/processed/alignments/qualimap/Mut4.part_007.vs.WT.part_001/qualimapReport.html"

# for i in {1..9}; do for file in $(ls OkinawaSeq/data/processed/vcf/hifiasm/Mut1.part_00*.vs.WT.part_00$i.ann.vcf.summary.genes.txt); do cmp $file OkinawaSeq/data/processed/vcf/hifiasm/Mut1.part_001.vs.WT.part_00$i.ann.vcf.summary.genes.txt; done ; done
# rule to summarize the differences using cmp
# for i in {1..9}; do for file in $(ls OkinawaSeq/data/processed/vcf/hifiasm/Mut1.part_00*.vs.WT.part_00$i.ann.vcf.summary.genes.txt); do cmp $file OkinawaSeq/data/processed/vcf/hifiasm/Mut1.part_001.vs.WT.part_00$i.ann.vcf.summary.genes.txt; done ; done


rule mutation_id:
    input:
        "data/processed/vcf/hifiasm/{sample}.vs.{ref}.ann.vcf".format(sample="{sample}", ref="WT")
    output:
        "data/processed/vcf/hifiasm/{sample}.vs.{ref}.ann.vcf.mutid".format(sample="{sample}", ref="WT")
    shell:
        "awk '!/^#/ {print $0}' {input} | awk -F \"|\" '{print $4, $10}' > {output}"

rule all:
    input:
        "data/processed/vcf/hifiasm/Mut1.part_004.vs.WT.part_001.q20.m1000.ann.vcf.summary.html",
        "data/processed/vcf/hifiasm/Mut1.part_004.vs.WT.part_001.q30.m1000.ann.vcf.summary.html",
        "data/processed/vcf/hifiasm/Mut2.part_004.vs.WT.part_001.q10.m2000.ann.vcf.summary.html",

        expand("data/processed/alignments/{sample}.part_00{p1}.vs.{sample}.part_00{p1}.asm.aln.sorted.bam.depth.png",sample=samples,p1=[1,2]),
        # expand("data/processed/igv/hifiasm/{sample}.part_00{p1}.vs.{sample}.part_00{p2}.delta.igv.html", sample="WT", p1=[1], p2=[4]),
        expand("data/processed/vcf/hifiasm/Mut{m}.part_00{p1}.vs.WT.part_00{p2}.ann.vcf.summary.html", m=[1, 2, 3, 4, 5], p1=[1], p2=[1, 2]),
        expand("data/processed/igv/hifiasm/Mut{m}.part_00{p1}.vs.WT.part_00{p2}.igv.html",             m=[1, 2, 3, 4, 5], p1=[1], p2=[1, 2]),

        # "data/processed/igv/hifiasm/Mut1.part_004.vs.Mut1.part_001.igv.html",
        # expand("data/processed/alignments/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.asm.aln.sorted.bam.bai", sample=["Mut1"], p1=[1],ref=["WT"],p2=[1]),
        # expand("data/processed/vcf/hifiasm/Mut4.part_001.vs.WT.part_00{p2}.ann.vcf.summary.html", p2=[1, 3, 4, 6, 7, 8]),
        # expand("data/processed/vcf/hifiasm/Mut3.part_001.vs.WT.part_00{p2}.ann.vcf.summary.html", p2=[1, 3, 4, 6, 7, 8]),
        # expand("data/processed/vcf/hifiasm/Mut2.part_005.vs.WT.part_00{p2}.ann.vcf.summary.html", p2=[1, 3, 4, 6, 7, 8]),

        # expand("data/processed/vcf/hifiasm/Mut5.part_006.vs.WT.part_00{p2}.ann.vcf", p2=[1, 3, 4, 6, 7, 8]),
        # expand("data/processed/vcf/hifiasm/Mut5.part_005.vs.WT.part_00{p2}.ann.vcf", p2=[1, 3, 4, 6, 7, 8]),
        # expand("data/processed/vcf/hifiasm/Mut5.part_005.vs.WT.part_00{p2}.ann.vcf", p2=[1, 3, 4, 6, 7, 8]),
        # expand("data/processed/vcf/hifiasm/Mut5.part_005.vs.WT.part_00{p2}.ann.vcf", p2=[1, 3, 4, 6, 7, 8]),
        # expand("data/processed/vcf/hifiasm/Mut5.part_005.vs.WT.part_00{p2}.ann.vcf", p2=[1, 3, 4, 6, 7, 8]),
        # expand("data/processed/vcf/hifiasm/Mut5.part_005.vs.WT.part_00{p2}.ann.vcf", p2=[1, 3, 4, 6, 7, 8]),

        # expand("data/processed/vcf/hifiasm/Mut{m}.part_00{p1}.vs.WT.part_00{p2}.ann.vcf", p2=[1, 3, 4, 6, 7, 8])

        # "data/processed/vcf/hifiasm/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.csq.vcf"
        # expand("{fastqc_html}",fastqc_html=glob.glob("data/processed/fastqc/*.html")),

        # expand("data/raw/fastq/{sample}.filtered.fastq.split/{sample}.filtered.part_00{p}.fastq.stats", sample=samples,p=split_numbers),
        # expand("data/raw/fastq/{sample}.filtered.fastq.split/{sample}.filtered.part_00{p}.fastq.stats", sample=samples, p=split_numbers),
        # expand("data/raw/fastq/{sample}.filtered.fastq.stats", sample=samples),
        # expand("data/raw/fastq/{sample}.fastq.stats", sample=samples),
        # expand("data/raw/fastq/{sample}_fastqc.html", sample=samples),
        # expand("data/processed/hifiasm/{sample}.filtered/{sample}.filtered.part_00{p}.asm.fasta.stats", sample=samples, p=split_numbers),

        # expand("data/processed/hifiasm/{sample}/{sample}.part_00{p}.asm.fasta.stats", sample=samples, p=split_numbers),
        # expand("data/processed/unicycler/{sample}.part_00{p}/assembly.fasta.stats", sample=samples, p=split_numbers),


# expand("data/processed/racon/{sample}.filtered.part_00{p}.racon.fastq.stats", sample=samples, p=split_numbers),
        # expand("data/processed/alignments/{sample}.vs.{sample}.split{s}.aln.sorted.bam.depth.png",sample=["WT","Mut1","Mut2","Mut3","Mut4","Mut5"], s=split_numbers),

        # expand("data/processed/assembly_comparison/{sample}/{sample}.split{i}vs{j}.dnadiff.png",
        #     sample=["WT", "Mut1", "Mut2", "Mut3", "Mut4", "Mut5"], i=split_numbers, j=split_numbers),

        # expand("data/processed/vcf/Mut{i}.split1.vs.WT.split1.vcf", i=range(1,6))

        # expand("data/processed/alignments/{sample}.vs.{sample}.split{s}.aln.sorted.bam.bai", sample
        #        =["WT", "Mut1", "Mut2", "Mut3", "Mut4", "Mut5"], s=split_numbers),

        # expand("{fasta}.stats",fasta=fasta_files)

rule bamtofastq:
    input: "data/raw/bam/{sample}.bam"
    output: "data/raw/fastq/{sample}.fastq"
    log: "logs/bamtofastq/{sample}.log"
    threads: 8
    shell: "samtools fastq -@ {threads} {input} > {output}  2>{log}"

rule fastqc:
    input:
        "{sample}.fastq"
    output:
        report("{sample}_fastqc.html", category="Raw Reads", subcategory="Quality Control")
    log: "logs/fastqc/{sample}.fastqc.log"
    shell: 'fastqc {input} 2>{log}'

rule filter:
    input:
        "{sample}.fastq"
    output:
        "{sample}.q{params.quality}.m{params.m}.filtered.fastq"
    log: "logs/fastq/{sample}.filter.log"
    params:
        quality=20,
        m=1000
    threads: 12
    shell: 'seqkit seq -Q {params.quality} -M {params.m} -j {threads} {input} > {output} 2>{log}'

rule split_fasta:
    input:
        "data/raw/fastq/{sample}.fastq"
    output:
        multiext("data/raw/fastq/{sample}.fastq.split/{sample}.", *expand("part_00{p}.fastq", p=split_numbers))
    log: "logs/split_fasta/{sample}.split.log"
    shell: 'seqkit split --quiet -w 0 --by-part {n_splits} {input} 2>{log}'

rule rule_graph:
    input:
        "{file}"
    output:
        "{file}.rulegraph.png"
    shell:
        "snakemake --dag {file} | dot -Tpng > {file}.rulegraph.png"

rule hifiasm:
    input:
        fastq="data/raw/fastq/{sample}.fastq.split/{sample}.part_00{p}.fastq"
    output:
        multiext(
            "data/processed/hifiasm/{sample}/{sample}.part_00{p}.",
            "a_ctg.gfa",
            "a_ctg.lowQ.bed",
            "a_ctg.noseq.gfa",
            "p_ctg.gfa",
            "p_ctg.lowQ.bed",
            "p_ctg.noseq.gfa",
        ),
    log:
        "logs/hifiasm/{sample}.part_00{p}.log",
    params:
        output_prefix="data/processed/hifiasm/{sample}/{sample}.part_00{p}",
        extra="--primary -f0 -k31 -w 128 -D20",
    threads: 64
    # shell: "hifiasm -o {params.output_prefix} -t {threads} {input.fastq} {params.extra} 2>{log}"
    shell: "hifiasm -o {params.output_prefix} -t {threads} {input.fastq} {params.extra} 2>{log}"

rule gfa_to_fasta:
    input:
        gfa = "data/processed/hifiasm/{sample}/{sample}.part_00{p}.p_ctg.gfa"
    output: "data/processed/hifiasm/{sample}/{sample}.part_00{p}.asm.fasta"
    log: "logs/assembly/{sample}.part_00{p}.gfa_to_fasta.log"
    shell: """awk '/^S/{{print \">\"$2\"\\n\"$3}}' {input.gfa} | seqkit replace -p ">*^" -r "{wildcards.sample}_c{{nr}}_" > {output} 2>{log}"""

rule minimap2_align_reads:
    input:
        reads="data/raw/fastq/{sample}.fastq.split/{sample}.part_00{p1}.fastq",
        reference="data/processed/hifiasm/{ref}/{ref}.part_00{p2}.asm.fasta",
        reference_index="data/processed/hifiasm/{ref}/{ref}.part_00{p2}.asm.fasta.amb"
    output:
        "data/processed/alignments/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.asm.aln.sam"
    log: "logs/alignments/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.asm.aln.sam.log"
    threads: 16
    shell: 'minimap2 -t {threads} -ax map-pb {input[1]} {input[0]} -o {output} 2>{log}'

rule compare_assemblies:
    input:
        "data/processed/hifiasm/{sample}/{sample}.part_00{p2}.asm.fasta",
        "data/processed/hifiasm/{sample}/{sample}.part_00{p1}.asm.fasta",
    output:
        multiext("data/processed/dnadiff/hifiasm/{sample}/{sample}.part_00{p1}.vs.part_00{p2}.dnadiff",
            ".report",
            ".delta",
            ".1delta")
    params:
        prefix="data/processed/dnadiff/hifiasm/{sample}/{sample}.part_00{p1}.vs.part_00{p2}.dnadiff"
    log: "logs/assembly/{sample}.{p1}.vs.{p2}.dnadiff.log"
    shell: 'dnadiff {input} -p {params.prefix} 2>{log}'

rule delta_plot:
    input:
        "{file}.1delta"
    output:
        report("{file}.1delta.png",category="Assembly crossvalidation")
    params:
        prefix = "{file}"
    log: "logs/assembly/{file}.mummerplot.log"
    shell: """
        mummerplot --postscript {input} -p {params.prefix} 2>{log}
        gnuplot {params.prefix}.gp 2>{log}
        gs -sDEVICE=png256 -sOutputFile={output} -r200 -dBATCH -dNOPAUSE -dEPSCrop -c "<</Orientation 2>> setpagedevice" -q {params.prefix}.ps -c quit 2>{log}
        """


rule dnadiff_vcf:
    input:
        delta = "data/processed/dnadiff/hifiasm/{sample}/{sample}.part_00{p1}.vs.part_00{p2}.dnadiff.1delta",
        query = "data/processed/hifiasm/{sample}/{sample}.part_00{p1}.asm.fasta",
        reference = "data/processed/hifiasm/{sample}/{sample}.part_00{p2}.asm.fasta",
    output:
        vcf=report("data/processed/dnadiff/hifiasm/{sample}/{sample}.part_00{p1}.vs.part_00{p2}.dnadiff.delta.vcf", category="Mutations", subcategory="vcf")
    log: "logs/vcf/{sample}.part_00{p1}.vs.part_00{p2}.vcf.log"
    shell: """delta2vcf < {input.delta} > {output.vcf}.full
    head -n 20 {output.vcf}.full > {output.vcf}"""

rule vcf:
    input:
        query = "data/processed/hifiasm/{sample}/{sample}.part_00{p1}.asm.fasta",
        reference = "data/processed/hifiasm/{ref}/{ref}.part_00{p2}.asm.fasta",
    output:
        report("data/processed/vcf/hifiasm/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.vcf", category="Mutations",subcategory="vcf")
    log: "logs/vcf/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.vcf.log"
    threads: 16
    shell: 'minimap2 -t {threads} -cx asm5 --cs {input.reference} {input.query} | sort -k6,6 -k8,8n  | paftools.js call -l 5 -L 6 -q 0 -f {input.reference} - > {output} 2>{log}'

rule snpeff:
    input:
        vcf = "{query}.vs.{reference}.vcf",
        gbk = "data/processed/snpeff/{reference}/genes.gbk",
    output:
        ann=report("{query}.vs.{reference}.ann.vcf", category="Mutations",subcategory="vcf"),
        html=report("{query}.vs.{reference}.ann.vcf.summary.html", category="Mutations", subcategory="Variant Analysis"),
        txt=report("{query}.vs.{reference}.ann.vcf.summary.genes.txt",category="Mutations",subcategory="affected genes")
    log: "logs/vcf/{query}.vs.{reference}.snpeff.vcf.log"
    shell: 'java -jar snpEff.jar ann {wildcards.reference} {input.vcf} \
    -c data/resources/snpEff.config -dataDir ../processed/snpeff -s {output.html} > {output.ann} 2>{log}'

rule igv_report:
    input:
        vcf = "data/processed/vcf/hifiasm/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.ann.vcf",
        fasta = "data/processed/hifiasm/{ref}/{ref}.part_00{p2}.asm.fasta",
        bam = "data/processed/alignments/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.asm.aln.sorted.bam",
        bai= "data/processed/alignments/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.asm.aln.sorted.bam.bai",
    output: report("data/processed/igv/hifiasm/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.igv.html",category="Mutations",subcategory="Alignment View")
    log: "logs/igv/{sample}.part_00{p1}.vs.{ref}.part_00{p2}.igv.log"
    shell:"create_report {input.vcf} \
            {input.fasta} \
            --flanking 1000 \
            --tracks  {input.bam}\
            --output {output} 2>{log}"


rule igv_report_compare:
    input:
        vcf = "data/processed/dnadiff/hifiasm/{sample}/{sample}.part_00{p1}.vs.part_00{p2}.dnadiff.delta.vcf",
        fasta = "data/processed/hifiasm/{sample}/{sample}.part_00{p1}.asm.fasta",
        bam = "data/processed/alignments/{sample}.part_00{p1}.vs.{sample}.part_00{p2}.asm.aln.sorted.bam",
        bai= "data/processed/alignments/{sample}.part_00{p1}.vs.{sample}.part_00{p2}.asm.aln.sorted.bam.bai",
    output: report("data/processed/igv/hifiasm/{sample}.part_00{p1}.vs.{sample}.part_00{p2}.delta.igv.html",category="Mutations",subcategory="Alignment View")
    log: "logs/igv/{sample}.part_00{p1}.vs.{sample}.part_00{p2}.igv.log"
    shell:"create_report {input.vcf} \
            {input.fasta} \
            --flanking 3000 \
            --tracks  {input.bam}\
            --output {output} 2>{log}"


rule prokka:
    input:
        # gb = "data/resources/CP020543.1.gbk",
        fasta = "data/processed/hifiasm/{sample}/{sample}.part_00{p}.asm.fasta"
    output:
        multiext("data/processed/prokka/hifiasm/{sample}.part_00{p}/{sample}.part_00{p}",
            ".gff",".gbk")
    log: "logs/prokka/{sample}.part_00{p}.log"
    shell: "prokka --force --genus Escherichia --debug --outdir data/processed/prokka/hifiasm/{wildcards.sample}.part_00{wildcards.p} --prefix {wildcards.sample}.part_00{wildcards.p} {input} 2>{log}"

rule snpeff_build_db:
    input:
        reference = "data/processed/hifiasm/{ref}/{ref}.part_00{p2}.asm.fasta",
        gbk = "data/processed/prokka/hifiasm/{ref}.part_00{p2}/{ref}.part_00{p2}.gbk"
    output:
        gbk="data/processed/snpeff/{ref}.part_00{p2}/genes.gbk",
        fasta="data/processed/snpeff/{ref}.part_00{p2}/{ref}.part_00{p2}.fasta"
    log:
        "logs/snpeff/{ref}.part_00{p2}.log"
    shell:
        """cp {input.reference} {output.fasta} 2>{log}
        cp {input.gbk} {output.gbk} 2>{log}
        echo "{wildcards.ref}.part_00{wildcards.p2}.genome: {wildcards.ref}.part_00{wildcards.p2}\n" >> data/resources/snpEff.config 2>{log}
        java -jar snpEff.jar build -c data/resources/snpEff.config -dataDir ../processed/snpeff -genbank -v {wildcards.ref}.part_00{wildcards.p2} 2>{log}"""


rule sam_to_bam:
    input:
        "data/processed/alignments/{file}.sam"
    output:
        "data/processed/alignments/{file}.bam"
    log: "logs/alignments/{file}.bam.log"
    shell: 'samtools view -b {input} > {output} 2>{log}'

rule bam_sort:
    input:
        "{file}.bam"
    output:
        "{file}.sorted.bam"
    log: "logs/alignments/{file}.sorted.bam.log"
    shell: 'samtools sort {input} -o {output} 2>{log}'

rule bam_index:
    input:
        "{file}.sorted.bam"
    output:
        "{file}.sorted.bam.bai"
    log: "logs/alignments/{file}.bam.bai.log"
    shell: 'samtools index {input} 1> {log}'

rule bam_depth:
    input:
        "{file}.sorted.bam"
    output:
        depth="{file}.sorted.bam.depth"
    log: "logs/alignments/{file}.bam.depth.log"
    shell: """samtools depth {input} > {output.depth} 2>{log}"""

rule bam_depth_plot:
    input:
        "{file}.bam.depth"
    output:
        png=report("{file}.bam.depth.png",category="Raw Reads", subcategory = "Coverage")
    log: "logs/alignments/{file}.bam.depth.png.log"
    shell: "python scripts/depth.py {input} {output.png} 2>{log}"

rule fasta_index:
    input:
        "{file}.fasta"
    output:
        "{file}.fasta.bwt",
        "{file}.fasta.amb",
        "{file}.fasta.ann",
        "{file}.fasta.pac",
        "{file}.fasta.sa"
    log: "logs/bwa_index/{file}.fasta.log"
    shell: 'bwa index {input} 2>{log}'

rule fasta_first10k:
    input:
        "{file}.fastq"
    output:
        "{file}.head.fastq"
    log: "logs/seqkit/{file}.head.fasta.log"
    shell: 'seqkit head -n 10000 {input} > {output} 2>{log}'

rule seqkit_fasta_stats:
    input:
        "{fasta}"
    output:
        "{fasta}.stats"
    log: "logs/assembly_stats/{fasta}.stats.log"
    shell: 'seqkit stats {input} > {output} 2>{log}'

rule qualimap:
    input:
        bam = "{file}.sorted.bam",
        fasta = "{file}.fasta"
    output:
        report = report("{file}.qualimap.html",category="Raw Reads", subcategory = "Alignment Quality"),
        bamqc = "{file}.bamqc"
    log: "logs/qualimap/{file}.qualimap.log"
    shell: "qualimap bamqc -bam {input.bam} -outdir {output.bamqc} -nt 32 -gd {input.fasta} 2>{log}"


# rule bam_quality_stats:
#     input:
#         bam = "{file}.sorted.bam",
#         fasta = "{file}.fasta"
#     output:
#         report = report("{file}.bam_quality_stats.html",category="Raw Reads", subcategory = "Alignment Quality"),
#         bamqc = "{file}.bamqc"
#     log: "logs/bam_quality_stats/{file}.bam_quality_stats.log"
#     shell: "bamqc -bam {input.bam} -outdir {output.bamqc} -nt 32 -gd {input.fasta} 2>{log}"

# samtools stats data/processed/alignments/Mut1.part_001.vs.Mut1.part_001.asm.aln.sorted.bam > samtoolsstats.txt
rule samtools_stats_bam:
    input:
        bam = "{file}.bam"
    output:
        stats = "{file}.bam.stats.txt"
    log: "logs/samtools_stats_bam/{file}.samtools_stats_bam.log"
    shell: "samtools stats {input.bam} > {output.stats} 2>{log}"