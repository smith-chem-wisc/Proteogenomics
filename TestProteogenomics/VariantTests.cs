﻿using NUnit.Framework;
using Proteogenomics;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TestProteogenomics
{
    [TestFixture]
    public class VariantTests
    {
        [OneTimeSetUp]
        public void Setup()
        {
            (new GeneModelTests()).Setup();
        }

        //[Test]
        //public void OneTranscriptOneHomozygous()
        //{
        //    Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));
        //    VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", "chr_1_one_homozygous_missense.vcf"));
        //    List<Variant> variants = vcf.Select(x => new Variant(null, x, genome)).ToList();
        //    Assert.AreEqual(1, variants.Count);

        //    GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_one_transcript.gtf"));
        //    List<Protein> proteins_wo_variant = geneModel.Translate(true).ToList();
        //    List<Transcript> transcripts = geneModel.ApplyVariants(variants);
        //    List<Protein> proteins = transcripts.Select(t => t.Protein()).ToList();
        //    Assert.AreEqual(1, geneModel.Genes.Count);
        //    Assert.AreEqual(1, proteins.Count);
        //    Assert.AreEqual(1, proteins_wo_variant.Count);
        //    Assert.AreEqual(2, new HashSet<string> { proteins[0].BaseSequence, proteins_wo_variant[0].BaseSequence }.Count);
        //    Assert.IsTrue(proteins[0].FullName != null);
        //    Assert.IsTrue(proteins[0].FullName.Contains(FunctionalClass.MISSENSE.ToString())); // sav
        //    Assert.IsTrue(proteins[0].FullName.Contains(GenotypeType.HOMOZYGOUS_ALT.ToString())); // sav
        //    Assert.IsTrue(proteins[0].FullName.Contains("1:69640"));

        //    string proteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr_1_one_homozygous_missense.fasta");
        //    ProteinDbWriter.WriteFastaDatabase(proteins, proteinFasta, " ");
        //    string[] proteinFastaLines = File.ReadLines(proteinFasta).ToArray();
        //    Assert.IsTrue(proteinFastaLines[0].Contains(FunctionalClass.MISSENSE.ToString())); //sav
        //    Assert.IsTrue(proteinFastaLines[0].Contains("1:69640"));
        //}

        //[Test]
        //public void OneTranscriptOneHeterozygous()
        //{
        //    Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));
        //    VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", "chr_1_one_heterozygous_missense.vcf"));
        //    List<Variant> variants = vcf.Select(x => new Variant(null, x, genome)).ToList();
        //    Assert.AreEqual(1, variants.Count);

        //    GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_one_transcript.gtf"));
        //    List<Protein> proteins_wo_variant = geneModel.Translate(true).ToList();
        //    List<Transcript> transcripts = geneModel.ApplyVariants(variants);
        //    List<Protein> proteins = transcripts.Select(t => t.Protein()).ToList();
        //    Assert.AreEqual(1, geneModel.Genes.Count);
        //    Assert.AreEqual(2, proteins.Count);
        //    Assert.AreEqual(1, proteins_wo_variant.Count);
        //    Assert.AreEqual(2, new HashSet<string> { proteins[0].BaseSequence, proteins[1].BaseSequence, proteins_wo_variant[0].BaseSequence }.Count);
        //    Assert.IsTrue(proteins.All(p => p.FullName != null));
        //    Assert.IsTrue(proteins.Any(p => p.FullName.Contains(FunctionalClass.MISSENSE.ToString()))); // sav
        //    Assert.IsTrue(proteins.Any(p => p.FullName.Contains(GenotypeType.HETEROZYGOUS.ToString()))); // sav
        //    Assert.IsTrue(proteins.Any(p => p.FullName.Contains("1:69640")));

        //    string proteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr_1_one_heterozygous_missense.fasta");
        //    ProteinDbWriter.WriteFastaDatabase(proteins, proteinFasta, " ");
        //    string[] proteinFastaLines = File.ReadLines(proteinFasta).ToArray();
        //    Assert.IsTrue(proteinFastaLines.Any(x => x.Contains(FunctionalClass.MISSENSE.ToString()))); // sav
        //    Assert.IsTrue(proteinFastaLines.Any(x => x.Contains("1:69640")));
        //}

        //[Test]
        //public void OneTranscriptOneHeterozygousSynonymous()
        //{
        //    Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));
        //    VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", "chr_1_one_heterozygous_synonymous.vcf"));
        //    List<Variant> variants = vcf.Select(x => new Variant(null, x, genome)).ToList();
        //    Assert.AreEqual(1, variants.Count);

        //    GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_one_transcript.gtf"));
        //    List<Protein> proteins_wo_variant = geneModel.Translate(true).ToList();
        //    List<Transcript> transcripts = geneModel.ApplyVariants(variants);
        //    List<Protein> proteins = transcripts.Select(t => t.Protein()).ToList();
        //    Assert.AreEqual(1, geneModel.Genes.Count);
        //    Assert.AreEqual(1, proteins.Count);
        //    Assert.AreEqual(1, proteins_wo_variant.Count);
        //    Assert.AreEqual(1, new HashSet<string> { proteins[0].BaseSequence, proteins_wo_variant[0].BaseSequence }.Count);
        //    Assert.IsTrue(proteins.Any(p => p.FullName.Contains(FunctionalClass.SILENT.ToString()))); // synonymous
        //    Assert.IsTrue(proteins.Any(p => p.FullName.Contains(GenotypeType.HETEROZYGOUS.ToString()))); // synonymous
        //    Assert.IsTrue(proteins.Any(p => p.FullName.Contains("1:69666")));

        //    string proteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", "chr_1_one_heterozygous_synonymous.fasta");
        //    ProteinDbWriter.WriteFastaDatabase(proteins, proteinFasta, " ");
        //    string[] proteinFastaLines = File.ReadLines(proteinFasta).ToArray();
        //    Assert.IsTrue(proteinFastaLines[0].Contains(FunctionalClass.SILENT.ToString())); // synonymous
        //    Assert.IsTrue(proteinFastaLines[0].Contains("1:69666"));
        //}

        //[Test]
        //public void OneTranscriptOneHomozygousSynonymous()
        //{
        //    Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));
        //    VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", "chr_1_one_homozygous_synonymous.vcf"));
        //    List<Variant> variants = vcf.Select(x => new Variant(null, x, genome)).ToList();
        //    Assert.AreEqual(1, variants.Count);

        //    GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_one_transcript.gtf"));
        //    List<Protein> proteins_wo_variant = geneModel.Translate(true).ToList();
        //    List<Transcript> transcripts = geneModel.ApplyVariants(variants);
        //    List<Protein> proteins = transcripts.Select(t => t.Protein()).ToList();
        //    Assert.AreEqual(1, geneModel.Genes.Count);
        //    Assert.AreEqual(1, proteins.Count);
        //    Assert.AreEqual(1, proteins_wo_variant.Count);
        //    Assert.AreEqual(1, new HashSet<string> { proteins[0].BaseSequence, proteins_wo_variant[0].BaseSequence }.Count);
        //    Assert.IsTrue(proteins.All(p => p.FullName != null));
        //    Assert.IsTrue(proteins.Any(p => p.FullName.Contains(FunctionalClass.SILENT.ToString()))); // synonymous
        //    Assert.IsTrue(proteins.Any(p => p.FullName.Contains(GenotypeType.HOMOZYGOUS_ALT.ToString()))); // synonymous
        //    Assert.IsTrue(proteins.Any(p => p.FullName.Contains("1:69666")));

        //    string proteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr_1_one_homozygous_synonymous.fasta");
        //    ProteinDbWriter.WriteFastaDatabase(proteins, proteinFasta, " ");
        //    string[] proteinFastaLines = File.ReadLines(proteinFasta).ToArray();
        //    Assert.IsTrue(proteinFastaLines[0].Contains(FunctionalClass.SILENT.ToString())); // synonymous
        //    Assert.IsTrue(proteinFastaLines[0].Contains("1:69666"));
        //}

        //[Test]
        //public void ProblematicChr19Gene()
        //{

        //    Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.19.fa"));
        //    GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "ProblematicChr19", "problematicChr19Gene.gff3"));
        //    geneModel.ApplyVariants(new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "ProblematicChr19", "chr19problematic.vcf")).Select(v => new Variant(null, v, genome.Chromosomes[0])).ToList());
        //}


        //[Test]
        //public void Chr19VariantTranscript()
        //{
        //    Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.19.fa"));
        //    GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "ProblematicChr19", "chr19variantTranscript.gff3"));
        //    var variants = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "ProblematicChr19", "chr19problematic.vcf"))
        //        .Select(v => new Variant(null, v, genome.Chromosomes[0]))
        //        .Where(v => v.SecondAlleleString.Length == 1 && v.ReferenceAlleleString.Length == 1).ToList();
        //    List <Transcript> transcripts = geneModel.ApplyVariants(variants).ToList();
        //    List<Protein> proteins = transcripts.Select(t => t.Protein(null)).ToList();
        //}

        // test todo: transcript with zero CodingSequenceExons and try to translate them to check that it doesn fail
        // test todo: multiple transcripts
    }
}