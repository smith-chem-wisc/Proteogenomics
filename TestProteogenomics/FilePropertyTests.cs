using NUnit.Framework;
using Proteogenomics;
using System.IO;
using System;

namespace TestProteogenomics
{
    [TestFixture]
    public class FilePropertyTests
    {
        [OneTimeSetUp]
        public void Setup()
        {
            (new GeneModelTests()).Setup();
        }

        [Test]
        public void CountReads()
        {
            var fastq = new FastqProperties(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper.fastq"));
            Assert.AreEqual(3970, fastq.ReadCount);
            Assert.AreEqual(65.3, Math.Round(fastq.AverageReadLength, 1));
        }

        [Test]
        public void BAMPropertiesStrandSpecificityTest()
        {
            BAMProperties bam = new BAMProperties(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "unstrandedSingle202122.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa")),
                0.8);
            Assert.AreEqual(Strandedness.None, bam.Strandedness);
            Assert.AreEqual(RnaSeqProtocol.SingleEnd, bam.Protocol);
        }
    }
}