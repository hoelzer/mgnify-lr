process sourmash_download_db {
        if (params.cloudProcess) { publishDir "${params.cloudDatabase}/sourmash", mode: 'copy', pattern: "gtdb.lca.json" }
        else { storeDir 'nextflow-autodownload-databases/sourmash' }  
	      label 'basics'
      output:
        file("gtdb.lca.json")
      script:
        """
        wget https://ndownloader.figshare.com/files/18809423?private_link=ed98a281ef089c033352 -O gtdb.lca.json
        """
    }

/* Comments:
This is the "autodownload" process for the sourmash databases.
It features auto storage via cloud (line 2) or local storage (line 3)

We use "storeDir" in local mode - because its like a "persistant" process output,
this is automatically always checked by nextflow. Its not usable for the cloud.
We use the "checkifexists" instead. (see main.nf code)
*/