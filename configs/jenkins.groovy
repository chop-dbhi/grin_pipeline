node ('respublica-slave') {
 
    //def workspace = "/data/svc_cbmicid/.jenkins/workspace/grin_test"
    def workspace = "/home/svc_cbmicid/.jenkins/workspace/grin_master"
    def workdir = "/mnt/isilon/cbmi/variome/zhangs3/projects/data/jenkins"
    def repodir = "/mnt/isilon/cbmi/variome/zhangs3/projects/repo/jenkins"
    env.PATH = "${workspace}/miniconda/bin:${env.PATH}"
    env.PATH = "${workspace}/miniconda/envs/grinenv/bin:${env.PATH}"
    
    stage 'Setup'
    sh """
        # cd ${workspace}
        # hostname
        # pwd
        #chmod og+rx /home/svc_cbmicid/.jenkins
        #chmod og+rx /home/svc_cbmicid/.jenkins/workspace
        #chmod -R og+r ${workspace}
        #rm -rf miniconda
        #./install_miniconda.sh ${workspace}
        #./setup_conda_env.sh ${workspace}
    """
    
    stage 'Cleanup'
    sh """
        #pwd
        #ls
        # env
        # which snakemake
        # which python
        # python --version
        cd ${workdir}
        #pwd
        #rm -rf lustre/*
        rm -rf lustre/GRCh37/picard
        rm -rf lustre/GRCh38/picard
        rm -rf GRCh37
        rm -rf GRCh38
        rm -rf .snakemake
    """

    stage 'Checkout' // need only do once
    // git url: 'https://github.research.chop.edu/BiG/grin.git' , branch: 'master'
    // def pullCommand = ["git", "pull", "origin", "master"]
    // Process pullProcess2 = pullCommand.execute() 
    sh """
        # source activate grinenv
        # pwd
        # ls
        cd ${repodir}
        git checkout master
        git pull
        chmod go+w Snakefile
        # cp ${workspace}/ops/*.* ${workspace}
        # chmod 755 ${workspace}/*.sh
        # rebuild environent if requirements file is changed
        diff -qB requirements.txt requirements.prev || conda install --name grinenv --file requirements.txt && cp requirements.txt requirements.prev
        # cp novoalign.lic ${workspace}/miniconda/envs/grinenv/bin
     """
    
    stage 'Build'
    // main point is to start snakemake in the working directory
    // env.PATH = "${workspace}/miniconda/bin:${env.PATH}"
    // env.PATH = "${workspace}/miniconda/envs/grinenv/bin:${env.PATH}"
    sh """
        cd ${workspace}
        # ls -R /home/svc_cbmicid/.jenkins/workspace/grin_test/ops
        # cp /home/svc_cbmicid/.jenkins/workspace/grin_test/ops/* .
        # cp /home/svc_cbmicid/.jenkins/workspace/grin_test/novoalign.lic .
        # cp ${workspace}/ops/*.* ${workspace}
        # chmod 755 ${workspace}/*.sh
        # ./install_miniconda.sh ${workspace}
        # cp /home/svc_cbmicid/.jenkins/workspace/grin_test/requirements.txt .
        # ls -lR miniconda
        # ./setup_conda_env.sh ${workspace}
        # ls -l miniconda
        # cp novoalign.lic ${workspace}/miniconda/envs/grinenv/bin
        source activate grinenv
        #which novoalign
        #which conda
        # env | sort
        # ls -R /home/svc_cbmicid/.jenkins/workspace/grin_test/miniconda/envs/grinenv
        #ls -l /home/svc_cbmicid/.jenkins/workspace/

        #ls /home/svc_cbmicid/.jenkins/workspace/grin_test/miniconda/envs
        #ls /home/svc_cbmicid/.jenkins/workspace/grin_test/miniconda/envs/grinenv/bin
        
        cd ${workdir}
        # whoami
        # hostname
        # who
        # hostname
        # which perl
        # perl -v
        snakemake --latency-wait 10 -j 10 --cluster-config configs/cluster.yaml -c 'qsub -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -l m_mem_free={cluster.m_mem_free} -pe smp {threads}' comvcfs
        # snakemake -j 10 comvcfs
        # snakemake dummy
        # ls -l /home/svc_cbmicid
        # cp /home/svc_cbmicid/snakejob.run_vep.0.sh.e4175621 .
        # file /mnt/isilon/cbmi/variome/vep/homo_sapiens/84_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.index
        # cp  /mnt/isilon/cbmi/variome/vep/homo_sapiens/84_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.index .
    """
}



