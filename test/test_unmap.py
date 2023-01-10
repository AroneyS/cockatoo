#!/usr/bin/env python3

import unittest
import os
import gzip
from bird_tool_utils import in_tempdir
import extern

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_conda = os.path.join(path_to_data,'.conda')

SAMPLE_READS_FORWARD = " ".join([
    os.path.join(path_to_data, "sample_1.1.fq"),
    os.path.join(path_to_data, "sample_2.1.fq"),
    os.path.join(path_to_data, "sample_3.1.fq"),
])
SAMPLE_READS_REVERSE = " ".join([
    os.path.join(path_to_data, "sample_1.2.fq"),
    os.path.join(path_to_data, "sample_2.2.fq"),
    os.path.join(path_to_data, "sample_3.2.fq"),
])
GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2.fna"),
    ])

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
APPRAISE_BINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "binned.otu_table.tsv")
APPRAISE_UNBINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "unbinned.otu_table.tsv")
ELUSIVE_CLUSTERS = os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters.tsv")


class Tests(unittest.TestCase):
    def test_unmap(self):
        with in_tempdir():
            cmd = (
                f"cockatoo unmap "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            bins_reference_path = os.path.join("test", "coassemble", "mapping", "sample_1_reference.fna")
            self.assertFalse(os.path.exists(bins_reference_path))

            output_bam_files = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.bam")
            self.assertFalse(os.path.exists(output_bam_files))

            coverm_working_dir = os.path.join("test", "coassemble", "mapping", "sample_1_coverm")
            self.assertFalse(os.path.exists(coverm_working_dir))

            unmapped_sample_1_path = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")
            self.assertTrue(os.path.exists(unmapped_sample_1_path))
            with gzip.open(unmapped_sample_1_path) as f:
                file = f.read().decode()
                self.assertTrue("@A00178:112:HMNM5DSXX:4:1622:16405:19194" in file)
                self.assertTrue("@A00178:118:HTHTVDSXX:1:1249:16740:14105" not in file)

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            test_dir = os.path.abspath("test")
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble -1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble"),
                        "-n 16 -t 16 -m 250 &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_assemble.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(coassemble_path) as f:
                self.assertEqual(expected, f.read())

            recover_path = os.path.join("test", "coassemble", "commands", "recover_commands.sh")
            self.assertTrue(os.path.exists(recover_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary recover --assembly", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta"),
                        "-1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_3_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_3_unmapped.2.fq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "recover"),
                        "-n 16 -t 16 -m 250 &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_recover.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

    def test_unmap_specified_files(self):
        with in_tempdir():
            cmd = (
                f"cockatoo unmap "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--appraise-binned {APPRAISE_BINNED} "
                f"--appraise-unbinned {APPRAISE_UNBINNED} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" not in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("target_elusive" not in output)
            self.assertTrue("cluster_graph" not in output)
            self.assertTrue("collect_genomes" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("aviary_commands" in output)

    def test_unmap_sra_download(self):
        with in_tempdir():
            cmd = (
                f"cockatoo unmap "
                f"--forward SRR8365216 SRR8365218 "
                f"--sra "
                f"--genomes {GENOMES} "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output_comb = extern.run(cmd)

            output_sra = output_comb.split("Building DAG of jobs...")[1]
            self.assertTrue("download_sra" in output_sra)
            self.assertTrue("aviary_commands" not in output_sra)

            output = output_comb.split("Building DAG of jobs...")[2]
            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" not in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("target_elusive" not in output)
            self.assertTrue("cluster_graph" not in output)
            self.assertTrue("collect_genomes" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("aviary_commands" in output)

    @unittest.skip("Downloads SRA data using Kingfisher. Test manually")
    def test_unmap_sra_download_real(self):
        with in_tempdir():
            cmd = (
                f"cockatoo unmap "
                f"--forward SRR9974654 SRR9974657 "
                f"--sra "
                f"--genomes {GENOMES} "
                f"--appraise-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
                f"--appraise-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
                f"--elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra.tsv')} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            sra_1_path = os.path.join("test", "coassemble", "sra", "SRR9974654_1.fastq.gz")
            self.assertTrue(os.path.exists(sra_1_path))
            with gzip.open(sra_1_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SRR9974654.1 1/1" in file)
                self.assertTrue("@SRR9974654.1 1/2" not in file)

            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR9974654_2.fastq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR9974657_1.fastq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR9974657_2.fastq.gz")))


if __name__ == '__main__':
    unittest.main()
