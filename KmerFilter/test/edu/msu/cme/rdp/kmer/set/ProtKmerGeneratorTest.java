/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.kmer.set;

import edu.msu.cme.rdp.kmer.Kmer;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class ProtKmerGeneratorTest {

    @Test
    public void testGenerator() {
        String seq = "tmamrqcalygkggigkstttqnlvaalaemgkkvmivgcdpkadstrlilhskaqgtvmemaasagsvedleledvlqigfggvkcvesggpepgvgcagrgvitainfleeegaysddldfvfydvlgdvvcg";
        ProtKmerGenerator kmer = new ProtKmerGenerator(seq, 10);
        long[] expected = new long[]{
            536708426961142L,
            286171060117187L,
            150274669009005L,
            305189780917667L,
            758873734624355L,
            640061464284267L,
            215768533929315L,
            149193644682349L,
            270597002464692L,
            777804730971791L,
            119953440559599L,
            460810377379311L,
            109233287183858L,
            117765469355593L,
            390795298851111L,
            120550587966705L,
            479919094406692L,
            720712232060036L,
            544793289068679L,
            544886647558372L,
            547874119228545L,
            643473212674088L,
            324944482403587L,
            265124275331181L,
            602677462699437L,
            145380390057393L,
            148572854466088L,
            250731715544331L,
            142115549520241L,
            44097957277219L,
            285234726028390L,
            120311978167488L,
            472283580831763L,
            476375797662317L,
            607326736240036L,
            294157143356544L,
            405829332668436L,
            601639670121103L,
            112171027550709L,
            211773161094823L,
            21341713978603L,
            682934847315303L,
            461816884079850L,
            141441501601108L,
            22528423864973L,
            720909563679140L,
            551107900880018L,
            746954225521219L,
            258637172983919L,
            395090187587057L,
            257987027516968L,
            374285532644609L,
            718137976201256L,
            462417101587716L,
            160648461852804L,
            637151151919252L,
            122638538248836L,
            546733503434883L,
            606973507276916L,
            282853816536721L,
            44122874434081L,
            286032075047968L,
            145827146793991L,
            162869070037217L,
            708210613820455L,
            144741505402081L,
            128128545496096L,
            722413735347217L,
            599241394258471L,
            35426199946482L,
            7738491444811L,
            247631726233955L,
            42915891588194L,
            247408623979587L,
            35776619448419L,
            18951915506801L,
            606461296217645L,
            266463062640038L,
            645518656582865L,
            390398687484449L,
            107859024233524L,
            73789054944899L,
            109449944551523L,
            124698505120883L,
            612652443340385L,
            464579770567731L,
            229853869213283L,
            599924373769329L,
            57281544293923L,
            707109510562918L,
            109506201160900L,
            126498716620931L,
            670259211341941L,
            56196532932259L,
            672389146989681L,
            124354473659947L,
            601643436590447L,
            112291554569700L,
            215630025702539L,
            144761381425513L,
            128764578245922L,
            742766783341639L,
            124639023237345L,
            610749023067169L,
            403670321824801L,
            532551323124771L,
            153143737353316L,
            396999967935638L,
            319099998671572L,
            78100795906688L,
            247425655328768L,
            36321622622215L,
            36392017068256L,
            38644639341570L,
            110728552087633L,
            165613946276386L,
            796046653473878L,
            703694960626368L,
            240603191313L,
            7699302122023L,
            246377667904739L,
            2786025053280L,
            89152801704977L,
            601089840874033L,
            94576491644454L,
            774647918937283L
        };

        String[] expectedKmers = new String[]{
            "tmamrqcaly",
            "mamrqcalyg",
            "amrqcalygk",
            "mrqcalygkg",
            "rqcalygkgg",
            "qcalygkggi",
            "calygkggig",
            "alygkggigk",
            "lygkggigks",
            "ygkggigkst",
            "gkggigkstt",
            "kggigksttt",
            "ggigkstttq",
            "gigkstttqn",
            "igkstttqnl",
            "gkstttqnlv",
            "kstttqnlva",
            "stttqnlvaa",
            "tttqnlvaal",
            "ttqnlvaala",
            "tqnlvaalae",
            "qnlvaalaem",
            "nlvaalaemg",
            "lvaalaemgk",
            "vaalaemgkk",
            "aalaemgkkv",
            "alaemgkkvm",
            "laemgkkvmi",
            "aemgkkvmiv",
            "emgkkvmivg",
            "mgkkvmivgc",
            "gkkvmivgcd",
            "kkvmivgcdp",
            "kvmivgcdpk",
            "vmivgcdpka",
            "mivgcdpkad",
            "ivgcdpkads",
            "vgcdpkadst",
            "gcdpkadstr",
            "cdpkadstrl",
            "dpkadstrli",
            "pkadstrlil",
            "kadstrlilh",
            "adstrlilhs",
            "dstrlilhsk",
            "strlilhska",
            "trlilhskaq",
            "rlilhskaqg",
            "lilhskaqgt",
            "ilhskaqgtv",
            "lhskaqgtvm",
            "hskaqgtvme",
            "skaqgtvmem",
            "kaqgtvmema",
            "aqgtvmemaa",
            "qgtvmemaas",
            "gtvmemaasa",
            "tvmemaasag",
            "vmemaasags",
            "memaasagsv",
            "emaasagsve",
            "maasagsved",
            "aasagsvedl",
            "asagsvedle",
            "sagsvedlel",
            "agsvedlele",
            "gsvedleled",
            "svedleledv",
            "vedleledvl",
            "edleledvlq",
            "dleledvlqi",
            "leledvlqig",
            "eledvlqigf",
            "ledvlqigfg",
            "edvlqigfgg",
            "dvlqigfggv",
            "vlqigfggvk",
            "lqigfggvkc",
            "qigfggvkcv",
            "igfggvkcve",
            "gfggvkcves",
            "fggvkcvesg",
            "ggvkcvesgg",
            "gvkcvesggp",
            "vkcvesggpe",
            "kcvesggpep",
            "cvesggpepg",
            "vesggpepgv",
            "esggpepgvg",
            "sggpepgvgc",
            "ggpepgvgca",
            "gpepgvgcag",
            "pepgvgcagr",
            "epgvgcagrg",
            "pgvgcagrgv",
            "gvgcagrgvi",
            "vgcagrgvit",
            "gcagrgvita",
            "cagrgvitai",
            "agrgvitain",
            "grgvitainf",
            "rgvitainfl",
            "gvitainfle",
            "vitainflee",
            "itainfleee",
            "tainfleeeg",
            "ainfleeega",
            "infleeegay",
            "nfleeegays",
            "fleeegaysd",
            "leeegaysdd",
            "eeegaysddl",
            "eegaysddld",
            "egaysddldf",
            "gaysddldfv",
            "aysddldfvf",
            "ysddldfvfy",
            "sddldfvfyd",
            "ddldfvfydv",
            "dldfvfydvl",
            "ldfvfydvlg",
            "dfvfydvlgd",
            "fvfydvlgdv",
            "vfydvlgdvv",
            "fydvlgdvvc",
            "ydvlgdvvcg"
        };

        int index = 0;

        while (kmer.hasNext()) {
            Kmer temp = kmer.next();
            long val = temp.getPart(0);
            assertEquals(expectedKmers[index], temp.toString());
            assertEquals(expected[index], val);
            assertEquals(index + 1, kmer.getPosition());
            index++;
        }
    }

    @Test
    public void testInvalidKmer() {
        String str = "acgto";
        try {
            new ProtKmerGenerator(str, 25);
            fail("Kmer size of 25 should have triggered an exception");
        } catch (IllegalArgumentException e) {
        }

        try {
            new ProtKmerGenerator(str, 6);
            fail("Smaller sequence than kmer should have triggered an exception");
        } catch (IllegalArgumentException e) {
        }

        try {
            ProtKmerGenerator kmer = new ProtKmerGenerator(str, 2);
            assertEquals(134L, (long) kmer.next().getPart(0));
            assertEquals(1, kmer.getPosition());
            assertEquals(195L, (long) kmer.next().getPart(0));
            assertEquals(2, kmer.getPosition());
            assertEquals(111L, (long) kmer.next().getPart(0));
            assertEquals(3, kmer.getPosition());
            long l = kmer.next().getPart(0);
            fail("Next should've triggered an exception for invalid character");
        } catch (IllegalArgumentException e) {
        }
    }
}
