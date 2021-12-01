process checkSampleSheet {
    label "artic"
    cpus 1
    input:
        file "sample_sheet.txt"
    output:
        file "samples.txt"
    """
    check_sample_sheet.py sample_sheet.txt samples.txt
    """
}


/**
 * Load a sample sheet into a Nextflow channel to map barcodes
 * to sample names.
 *
 * @param samples CSV file according to MinKNOW sample sheet specification
 * @return A Nextflow Channel of tuples (barcode, sample name)
 */ 
def check_sample_sheet(samples)
{
    println("Checking sample sheet.")
    sample_sheet = Channel.fromPath(samples, checkIfExists: true)
    sample_sheet = checkSampleSheet(sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.barcode, row.sample_name) }
    return sample_sheet
}


/**
 * Find fastq data using various globs. Wrapper around Nextflow `file`
 * method.
 *
 * @param patten glob pattern for top level input folder.
 * @param maxdepth maximum depth to traverse
 * @return list of files.
 */ 
def find_fastq(pattern, maxdepth)
{
    files = []
    extensions = ["fastq", "fastq.gz", "fq", "fq.gz"]
    for (ext in extensions) {
        files += file("${pattern}/*.${ext}", type: 'file', maxdepth: maxdepth)
    }
    return files
}


/**
 * Rework EPI2ME flattened directory structure into standard form
 * files are matched on barcode\d+ and moved into corresponding
 * subdirectories ready for processing.
 *
 * @param input_folder Top-level input directory.
 * @param output_folder Top-level output_directory.
 * @return A File object representating the staging directory created
 *     under output_folder
 */ 
def sanitize_fastq(input_folder, output_folder)
{
    println("Running sanitization.")
    println(" - Moving files: ${input_folder} -> ${output_folder}")
    staging = new File(output_folder)
    staging.mkdirs()
    files = find_fastq("${input_folder}/**/", 1)
    for (fastq in files) {
        fname = fastq.getFileName()
        // find barcode
        pattern = ~/barcode\d+/
        matcher = fname =~ pattern
        if (!matcher.find()) {
            // not barcoded - leave alone
            fastq.renameTo("${staging}/${fname}")
        } else {
            bc_dir = new File("${staging}/${matcher[0]}")
            bc_dir.mkdirs()
            fastq.renameTo("${staging}/${matcher[0]}/${fname}")
        }
    }
    println(" - Finished sanitization.")
    return staging
}


/**
 * Resolves input folder containing barcode subdirectories
 * or a flat set of fastq data to a Nextflow Channel. Removes barcode
 * directories with no fastq files.
 *
 * @param input_folder Top level input folder to locate fastq data
 * @param sample_sheet List of tuples mapping barcode to sample name
 *     or a simple string for non-multiplexed data.
 * @param min_barcode minimum barcode of valid range (default: null)
 * @param max_barcode maximum barcode of valid range (default: null)
 * @return Channel of tuples (path, sample_name)
 */
def resolve_barcode_structure(input_folder, sample_sheet, min_barcode, max_barcode)
{
    println("Checking input directory structure.")
    barcode_dirs = file("$input_folder/barcode*", type: 'dir', maxdepth: 1)
    not_barcoded = find_fastq("$input_folder/", 1)
    samples = null
    //check min max are not null
    if (!min_barcode){
        min_barcode = 0
    }
    if (!max_barcode){
        max_barcode = Integer.MAX_VALUE
    }
    if (barcode_dirs) {
        println(" - Found barcode directories")
        // remove empty barcode_dirs
        valid_barcode_dirs = []
        invalid_barcode_dirs = []
        out_of_range_barcode_dir = []

        for (d in barcode_dirs) {
            // find barcode
            pattern = ~/barcode(\d+)/
            matcher = "${d}" =~ pattern
            value = matcher[0][1].toInteger()
            if(!find_fastq(d, 1)) {
                invalid_barcode_dirs << d
            }else if(value >= min_barcode && value <= max_barcode){
                valid_barcode_dirs << d 
            }else{
                out_of_range_barcode_dir << d
            }
         }
        
        if (invalid_barcode_dirs.size() > 0) {
            println(" - Some barcode directories did not contain .fastq(.gz) files:")
            for (d in invalid_barcode_dirs) {
                println("   - ${d}")
            }
        }
        if (out_of_range_barcode_dir.size() > 0) {
            println(" - Some barcode directories were outside of expected range")
            for (d in out_of_range_barcode_dir) {
                println("   - ${d}")
            }
        }
        // link sample names to barcode through sample sheet
        if (!sample_sheet) {
            sample_sheet = Channel
                .fromPath(valid_barcode_dirs)
                .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
                .map { path -> tuple(path.baseName, path.baseName) }
        }
        samples = Channel
            .fromPath(valid_barcode_dirs)
            .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
            .map { path -> tuple(path.baseName, path) }
            .join(sample_sheet)
            .map { barcode, path, sample -> tuple(path, sample) }
    } else if (not_barcoded) {
        println(" - Found fastq files, assuming single sample")
        sample = (sample_sheet == null) ? "unknown" : sample_sheet
        samples = Channel
            .fromPath(input_folder, type: 'dir', maxDepth:1)
            .map { path -> tuple(path, sample) }
    }
    return samples
}


/**
 * Take an input directory and sample sheet to return a channel of
 * named samples.
 *
 * @param input_folder Top level input folder to locate fastq data
 * @param sample_sheet List of tuples mapping barcode to sample name
 *     or a simple string for non-multiplexed data.
 * @param min_barcode minimum barcode of valid range (default: null)
 * @param max_barcode maximum barcode of valid range (default: null)
 * @return Channel of tuples (path, sample_name)
 */
def fastq_ingress(input_folder, output_folder, samples, sanitize, min_barcode, max_barcode)
{
    // EPI2ME harness 
    if (sanitize) {
        staging = "${output_folder}/staging"
        input_folder = sanitize_fastq(input_folder, staging)
    }
    // check sample sheet
    sample_sheet = null
    if (samples) {
        sample_sheet = check_sample_sheet(samples)
    }
    // resolve whether we have demultiplexed data or single sample
    if (!params.multiple_samples) {
        data = resolve_barcode_structure(input_folder, sample_sheet, min_barcode, max_barcode)
        // return error if data empty after processing
        if (data == null) {
            println("")
            println("Error: `--fastq` Unable to find FASTQ files or BARCODE folders in the provided --fastq path")
            exit 1
        }
    } else {
        pattern = "${input_folder}/*.fastq.gz"
        data = Channel.fromPath(pattern, checkIfExists: true)
                      .map { path -> tuple(path, path.simpleName) }
    }
    
    return data
}
