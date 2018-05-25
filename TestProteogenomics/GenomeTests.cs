using NUnit.Framework;
using Proteogenomics;
using System.IO;

namespace TestProteogenomics
{
    [TestFixture]
    public class GenomeTests
    {
        [Test]
        public void KaryotypicOrder()
        {
            Genome headers = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "headers.fa"));
            var seqs = headers.KaryotypicOrder();
            Assert.IsTrue(seqs[0].FriendlyName == "chr1" && seqs[1].FriendlyName == "chr2");
            Assert.IsTrue(seqs[22].FriendlyName == "chrX" && seqs[23].FriendlyName == "chrY" && seqs[24].FriendlyName == "chrM");
        }

        [Test]
        public void KaryotypicOrderShort()
        {
            Genome headers = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "headersShort.fa"));
            var seqs = headers.KaryotypicOrder();
            Assert.IsTrue(seqs[0].FriendlyName == "chr9" && seqs[1].FriendlyName == "chr20");
        }
    }
}