# citation_network_my_tools
 small scripts I used when develop citation network

操作流程：

第一步：复制网页上的文章数据到一个txt。

第二步：修改get_query_expression程序中的“文件路径”为你刚建立的txt的路径，然后运行程序。会生成一个后缀为.expr.txt的文件。

第三步：把.expr.txt文件的内容复制到web of science core collection的高级检索一栏，然后检索。

第四步：按照课上讲过的流程把检索出来这些文献的完整引文数据导出（此处见课堂ppt），重命名为download_*year*.txt。

第五步：按课上教的流程运行并设置citespace。注意这里右侧部分中间的节点类型选择，必须只选“article”，其他都不点。点击GO开始生成网络，第一个弹窗问是否按年份，我们只有一个年份，改为NO。其他直接默认就行。

第六步：生成网络后，不用任何操作，我们不用这软件分析。在上边栏选export->network->adjacent matrix，导出邻接矩阵。同样继续选择导出network summary。此时你的结果文件夹中会多出三个文件——
1.network*x*.txt 2.network*x*Labels.txt 3.summary.tsv
将summary文件重命名成“network*x*summary.tsv”便于后续操作。

第七步：修改PI_extraction程序中的文件路径为你存放download_*year*.txt文件的路径，并在路径下放入参考的PI清单（"pi_list.txt"），最后运行程序（注意运行时间较长，大约三十秒）。该路径下会生成名为pi-title.tsv的文件。确认将该文件移动到第六步三个文件同一路径下。

第八步：打开R代码，将开头的loc地址修改为你第六、七步四个文件的存储路径，将下面的网络名称修改为你的“network*x*”。然后ctrl+shift+enter运行全部。会在原目录下生成两个新文件：
1.network*x*_readable.tsv 2.network*x*_title.tsv

第九步：打开cytoscape，安装app：aMatloader（好像是叫这个，读取邻接矩阵的）和autoannotation（做聚类的）。用aMat读入邻接矩阵network*x*_readable.tsv，生成网络。这里可以根据edge的强度让edge的颜色渐变。接下来import table读入title注释network*x*_title.tsv，读的时候检查一下要把author那一栏标记为key

第十步：到这里网络应该已经生成好了，你可以把node label改成title。然后进行聚类的话比较简单，直接打开autoannotate，上半部分不用动，中间layout不重叠打勾，下面label选择title。然后运行，就会自动NSL聚类了
