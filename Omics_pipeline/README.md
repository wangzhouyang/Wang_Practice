# Omics_pipeline
This pipeline is built to ease my pressure for Multiple omics analysis. In this version, I'm focused on the process of data polishing.

# Get Start & Usage

本项目目前仅包含**Metagenomics宏基因组学标准流程**和**Lipidomics脂质组学分析流程**。

### Metagenomics pipeline

在集群工作目录下`clone`本仓库并将主程序添加到环境变量中。

```
cd /your/dir/
clone git@biogit.cn:Fangchao/Omics_pipeline.git
ln -s /your/dir/Omics_pipeline/MetaGenomics/cOMG ~/bin/
```

#### Usage:

```
cOMG 		#直接执行命令可以查看使用说明
usage:
        cOMG <pe|se> [options]
pattern
        pe|se           pair end | single end
options:
        -p|path         :[essential]sample path file (SampleID|fqID|fqPath)
        -i|ins          :[essential for pair-end seq]insert info file
        -s|step         :functions,default 1234
                             1       trim+filter, OA method
                             2       remove host genomic reads
                             3       soap mapping to microbiotic genomics
                             4       combine samples' abun into a single profile table
        -o|outdir       :output directory path. Conatins the results and scripts.
        -c|config       :provide a configure file including needed database and parameters for each setp
        -h|help         :show help info
        -v|version      :show version and author info.
```

**path file**: 用于记录raw data文件位置和id信息的文件，每行三列分别记录下**样本编号**, **数据编号** 和 **fq文件路径**。

- `样本编号`：生物学，统计学意义上的样本个体，用于后续分析的基本个体
- `数据编号`：如果同一个样本进行多次测序，则会产生多个数据，此时需要用数据编号来区分（可以使文库号，日期，批次，等等）。`注意`:拥有相同`样本编号`的多个数据会被最终合并计算相对丰度。
- `fastq路径`: 必须是工作环境可以访问到的路径位置

**config file**: 由于流程涉及到的分析步骤较多，对于每个具体工具的参数定义，统一放在配置文件中进行处理：

```
###### configuration

### Database location
db_host = $META_DB/human/hg19/Hg19.fa.index	#宿主参考基因集db前缀，用于去除宿主来源的reads
db_meta = $META_DB/1267sample_ICG_db/4Group_uniqGene.div_1.fa.index,$META_DB/1267sample_ICG_db/4Group_uniqGene.div_2.fa.index #参考基因集db前缀，多个文件可以用逗号分隔

### reference gene length file
RGL  = $META_DB/IGC.annotation/IGC_9.9M_update.fa.len #与参考基因集匹配的每个基因的长度信息，用于计算相对丰度
### pipeline parameters
PhQ = 33        		#reads Phred Quality system: 33 or 64.
mLen= 30                #minimal read length allowance
seedOA=0.9				
fragOA=0.8

qsub = 1234             #Following argment will enable only if qusb=on, otherwise you could commit it
q   = st.q              #queue for qsub
P   = st_ms             #Project id for qsub
B   = 1					#全局设定投递任务的备份数
B1  = 3					#针对第一步的任务投递备份数
p   = 6                 #全局计算核心数
p1  = 1					#具体到第一步的计算核心数，该参数比全局设定优先级高
p4  = 1					#具体到第四步的计算核心数，该参数比全局设定优先级高
f1  = 0.5G              #virtual free for qsub in step 1 (trim & filter)
f2  = 6G                #virtual free for qsub in step 2 (remove host genes)
f3  = 14G               #virtual free for qsub in step 3 (aligned to gene set)
f4  = 8G                #virtual free for qsub in step 4 (calculate soap results to abundance)
s   = 120				#qusbM定时检查任务完成情况的时间间隔（秒）
r   = 10                #repeat time when job failed or interrupted
```

上述配置文件准备完毕后，运行本脚本可以生成工作目录：

```
cOMG se -p sample.path.file -o demo
```

随后进入工作目录，检查脚本无误后可以启动执行脚本：

```
sh RUN.batch.sh			# 模式一，全部样本完成当前步骤后才会进入下一步骤；
sh RUN.linear.1234.sh	# 模式二，每个样本依次运行每个步骤，相互不影响；
sh RUN.qsubM.sh			# 模式三，同上，采用改进的qsub管理脚本，可自动处理异常情况
```

完成后可以执行`sh report.stat.sh`打印报告表格。

若中途出现错误，可以进入`script`目录对个别脚本进行调试。

--------------

### Lipidomics pipeline

### What do I wanna perform?
As a pipeline, I plan to orgainze the workshop directory like this:
```
./                 #output directory
|-- pip.work.sh		# A script contained all function set to run. A	MAIN SWITCH.
|-- Shell			# A directory contains all of scripts organized under steps order.
|-- Results			# A directory contains all of the results (exactly, including the intermediate results)
```
And all the users need to do is preparing a config file and write them into a script to build the workshop above.
Here is an example:

 `/ifs1/ST_MD/PMO/F13ZOOYJSY1389_T2D_MultiOmics/Users/fangchao/lipidomics.20151118/pip.config.sh`

For a better understanding of the pipeline's logic, a tree following shows you how the pip.work.sh works:
```
./pip.work.sh
	|--> sh step1.sh
			|--> sh function1.sh
					|--> sh/perl sub-function scripts/software/program [parameters]
			|--> sh function2.sh
					|--> sh/perl sub-function scripts/software/program [parameters]
			...
```
As you can see, the sub-funtion tools could come from websites, packages, or just wirtten by yourself. And what you need to do is to locate the scripts pathway and make sure the parameters are friendly for most of the naming manners, such as the capablility to read and locate an absolute path. Thus you can leave the rest things to the pipeline.

In the following step, I'll add your scripts into pipeline and distribute the unified input parameters as well as a proper output directory. Or some addtional options for the function of your part.


```

```