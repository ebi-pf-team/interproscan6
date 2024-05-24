# Understanding HMMER.OUT parser

<!-- TOC -->
- [Initial info](#initial-info)
  - [Query](#query)
  - [Accession](#accession)
  - [Description](#description)
- [List of hits](#list-of-hits)
- [Domain and Alignment Annotation](#domain-and-alignment-annotation)
  - [How we get the alignments](#how-we-get-the-alignments)
- [Two good examples](#two-good-examples)
<!-- /TOC -->

## [Initial info](initial-info)

### Query
Name of the query sequence that input sequence matched, followed by qlen.

Sample:

```Query:       exonuc_ExeM-GG  [M=884]```

```Query:       PTHR31901.orig.30.pir  [M=558]```


### Accession
Each hmmer.out file is generated for a member hmm, in this case it's the member accession that was matched.

Sample: 

```Accession:   NF033680.1```

```Accession:   ANF00001```

### Description
description/detailed name of the matched sequence

Sample:

```Description: extracellular exonuclease ExeM```

```Description: Translation of E. coli type CRISPR repeat```


## List of hits
List of the most significant sequence hits, ordered by E-value.

## Domain and Alignment Annotation
Details the location of the domains identified in each sequence.

#### How we get the alignments
On the I5 we had this regex to get the alignments:

```^\\s+(\\w+)\\s+(\\S+)\\s+([-a-zA-Z]+)\\s+(\\S+)\\s*$```

On I6 we change to:

```^\\s+(\\S+)\\s+(\\S+)\\s+([-a-zA-Z]+)\\s+(\\S+)\\s*$```

To be more specific, in I6 we have extended \w+ (restricted to letters, digits and underscores) to \S+ to not lose alignments that contain other characters types in the sequence id.
I5 deal with this by treating the fasta before sending it to the hmmer (changing temporarily the sequence id).

On the I5 parse this was necessary so that it could only get the alignment to the matched sequence, and not include the alignment (comparative) of the input sequence together.
But in the I6 parse, the way we deal with the current_sequence is slightly different, so we were able to have this regex most inclusive and just adding a check so that the id of the alignment besides having a match with the regex need to be equal to the current_sequence.



### Two good examples

To understand the hmmer.out matches and what the result of the parse looks like, we get two good examples:

- First one is a example with multiple domains:

```
Query:       PTHR46751.orig.30.pir  [M=266]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence Description
    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
    1.4e-10   52.7  30.0    8.5e-05   33.7  13.3    2.0  2  P22298


Domain annotation for each sequence (and alignments):
>> P22298
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   33.7  13.3   1.3e-12   8.5e-05      82     147 ..      23      96 ..       3      98 .. 0.81
   2 !   28.2  11.4     6e-11    0.0039      85     130 ..      80     126 ..      75     129 .] 0.90

  Alignments for each domain:
  == domain 1  score: 33.7 bits;  conditional E-value: 1.3e-12
  PTHR46751.orig.30.pir  82 lteslfprrCPkik.ekCefkerdeCtkdksCpekekCCvfsCGkkCldleq.......diCelPketGpClal 147
                              ++l  + CP  k ++C   e+ +Ct d +Cp+k+kCC   C+ kCl++            ++P   G C++l
                 P22298  23 AENALKGGACPPRKiVQCLRYEKPKCTSDWQCPDKKKCCRDTCAIKCLNPVAitnpvkvKPGKCPVVYGQCMML 96
                            55678889998777369*********************************998888888889999999999986 PP

  == domain 2  score: 28.2 bits;  conditional E-value: 6e-11
  PTHR46751.orig.30.pir  85 slfprrCPkikekCef.kerdeCtkdksCpekekCCvfsCGkkCldl 130
                            ++ p++CP +  +C + +  ++C++d +C    kCC+  CGk+Cl++
                 P22298  80 KVKPGKCPVVYGQCMMlNPPNHCKTDSQCLGDLKCCKSMCGKVCLTP 126
                            5789**********983568*************************98 PP
```
```json
{
  "P22298": {
    "PTHR46751": {
      "accession": "PTHR46751",
      "name": "PTHR46751.orig.30.pir",
      "description": "",
      "evalue": "1.4e-10",
      "score": "52.7",
      "qlen": "266",
      "bias": "30.0",
      "member_db": "panther",
      "version": "18.0",
      "model-ac": "PTHR46751",
      "locations": [
        {
          "start": "23",
          "end": "96",
          "representative": "",
          "hmmStart": "82",
          "hmmEnd": "147",
          "hmmLength": "266",
          "rawHmmBounds": "..",
          "hmmBounds": "Incomplete",
          "evalue": "8.5e-05",
          "score": "33.7",
          "envelopeStart": "3",
          "envelopeEnd": "98",
          "postProcessed": "false",
          "alignment": "AENALKGGACPPRKiVQCLRYEKPKCTSDWQCPDKKKCCRDTCAIKCLNPVAitnpvkvKPGKCPVVYGQCMML",
          "cigar_alignment": "14M1I37M7I15M"
        },
        {
          "start": "80",
          "end": "126",
          "representative": "",
          "hmmStart": "85",
          "hmmEnd": "130",
          "hmmLength": "266",
          "rawHmmBounds": "..",
          "hmmBounds": "Incomplete",
          "evalue": "0.0039",
          "score": "28.2",
          "envelopeStart": "75",
          "envelopeEnd": "129",
          "postProcessed": "false",
          "alignment": "KVKPGKCPVVYGQCMMlNPPNHCKTDSQCLGDLKCCKSMCGKVCLTP",
          "cigar_alignment": "16M1I30M"
        }
      ]
    }
  }
}
```

- On second one we have two queries (exonuc_ExeM-GG and ExeM_NucH_DNase) matching with two sequences (WP_249252303.1 and WP_338726824.1)

```
(...)
Query:       exonuc_ExeM-GG  [M=884]
Accession:   NF033680.1
Description: extracellular exonuclease ExeM
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence       Description
    ------- ------ -----    ------- ------ -----   ---- --  --------       -----------
          0 1701.4   8.9          0 1701.2   8.9    1.0  1  WP_249252303.1  extracellular exonuclease ExeM [Shewanella aq
          0 1609.3   7.1          0 1609.1   7.1    1.0  1  WP_338726824.1  extracellular exonuclease ExeM [Shewanella ba


Domain annotation for each sequence (and alignments):
>> WP_249252303.1  extracellular exonuclease ExeM [Shewanella aquimarina]
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ! 1701.2   8.9         0         0       1     884 []       1     885 []       1     885 [] 0.99

  Alignments for each domain:
  == domain 1  score: 1701.2 bits;  conditional E-value: 0
  exonuc_ExeM-GG   1 mdnvkklsvlalsiaaalpllasadvliseyveGssnnkaielynsgdaavdltGyklvrykdGataasdlvaldgqtlaaksvkvivhnsaait 95 
                     mdnvkklsvlalsiaaalp++a+ad+liseyveG+snnkaielynsgdaa+dltGyklvrykdGat+a d+valdg t++ak++kvivh+sa+i+
  WP_249252303.1   1 MDNVKKLSVLALSIAAALPVSAKADLLISEYVEGNSNNKAIELYNSGDAALDLTGYKLVRYKDGATNALDMVALDGLTIDAKKLKVIVHPSAEIA 95 
                     9********************************************************************************************** PP

  exonuc_ExeM-GG  96 ldagvdsltgnlyfnGddavallkdgvvvdvvGeiptpddwakdvtlqrkadvlaastvydesewttkakdtfsGlGarsgdfvggaatepevpa 190
                     l ++vds+tgnlyfnG+daval+kdg+vvd++G++ptp++w++++t+qr++d+l+as+v++es+w+ ++ dtfsGlG+++g   ++a+tepe+pa
  WP_249252303.1  96 LGSDVDSSTGNLYFNGGDAVALVKDGAVVDIIGAVPTPSGWGMNTTMQRNTDALSASNVFVESQWKVLPVDTFSGLGNLDG---SSAPTEPEAPA 187
                     ********************************************************************************7...56899****** PP

  exonuc_ExeM-GG 191 fscvGakltpiydvqGaGessplvaegkyesdeevtvkGvvtargdslfkGfylqevkGdnspytsdGifvylgeaapeaiqpGvevcvqGkvke 285
                     f+c+Gak+tpiy+vqGaG+sspl+a+g+y s+eevt++GvvtargdslfkGfylqev+Gd+spytsdGifvylgeaa e+iqpGvevcvqG vke
  WP_249252303.1 188 FTCTGAKITPIYSVQGAGKSSPLIADGEYASKEEVTIRGVVTARGDSLFKGFYLQEVQGDGSPYTSDGIFVYLGEAAGEEIQPGVEVCVQGLVKE 282
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 286 yygltqidikadkklevgakaevpaaealyvadGetladaleryeGmkvkldagsdmkvtrtfsydyasrrnnlvlshkaplmkptqvyaalsde 380
                     yygltq+dikadkk+evg++ ++paa+++yvadGe+l+daler+eGm+v+ldagsdmkv+rtfsydya+rrnn++lshkaplmkptq+y+a++de
  WP_249252303.1 283 YYGLTQLDIKADKKMEVGENLGAPAATPFYVADGENLGDALERFEGMNVVLDAGSDMKVSRTFSYDYAARRNNMMLSHKAPLMKPTQLYPAMTDE 377
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 381 aialekknkanqlfiesdlkaknGevpylpdfnaetGyirvGdqltnleGvvgysygayrlvatneitagdlirendrtdapevatkGdirvasf 475
                     aialekkn+anqlf+esd+kaknGevpy+p+fnaetGyirvGdqltnleGv++ysy++yrlvatnei+agd+ire+dr+dapevatkGdirvasf
  WP_249252303.1 378 AIALEKKNRANQLFVESDFKAKNGEVPYFPTFNAETGYIRVGDQLTNLEGVISYSYNEYRLVATNEIVAGDFIREDDRVDAPEVATKGDIRVASF 472
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 476 nvlnlftsddavgg......kddadadeasrgcgrGaltleeyllqrtkivnalkemnadivGlmeiennGfgeksaiqylldalnaelpaaday 564
                     nvlnlftsd++vgg      kd+adad asrgcgrGa+tle+yllqrtkiv+al+emnadi+GlmeiennGfge+s+iqyl+++lna+lpaaday
  WP_249252303.1 473 NVLNLFTSDSDVGGalnascKDQADAD-ASRGCGRGAHTLEDYLLQRTKIVSALTEMNADIIGLMEIENNGFGEDSSIQYLVNDLNANLPAADAY 566
                     *******************99999999.999**************************************************************** PP

  exonuc_ExeM-GG 565 kfieiadadkyegkfiGsdaitvGllyraakvkpegdafvietpeqhaaeGvatreeedkvetspaydkyqrhslaqtfkikdekltvvvnhlks 659
                     kfie+adadky+g+++G daitvG+lyr++kv+p+gdafvietpeqha+eGv+tr+e+dk+e spa+dkyqrhsl+q+f++++ekltvvvnhlks
  WP_249252303.1 567 KFIEVADADKYKGEYVGGDAITVGMLYRPSKVEPTGDAFVIETPEQHAMEGVVTRGEGDKLELSPAFDKYQRHSLGQAFSVNGEKLTVVVNHLKS 661
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 660 kGsgcledwvefeeskdpadlqGkcnefrvsaakvlGeslkdveGdllviGdlnayGledpvrvltdydaatserdivtaswttleGkvyereGs 754
                     kGsgc+edw++fe+s+dpadlqGkcnefrvsaakv+G++l+dveGdllviGd+nayGledp+rvltdydaa+s+rdi+taswttl GkvyereGs
  WP_249252303.1 662 KGSGCIEDWANFEDSSDPADLQGKCNEFRVSAAKVIGDALQDVEGDLLVIGDMNAYGLEDPIRVLTDYDAANSDRDIMTASWTTLGGKVYEREGS 756
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 755 kiekGyGlinlntqahGaatysysynGelGnldhalanaslaarlvdiedwhinsvesnlfeysskytGdleksenafsasdhdpvivalsypap 849
                     kiekGyG++nlntq+hG+atysysynGelGnldhalan+slaar+vdi+dwhinsvesnlfeysskytG+leksenafsasdhdpvivalsypap
  WP_249252303.1 757 KIEKGYGMVNLNTQVHGTATYSYSYNGELGNLDHALANTSLAARVVDIQDWHINSVESNLFEYSSKYTGELEKSENAFSASDHDPVIVALSYPAP 851
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 850 veevepepekkkddGGslGylalallsllGlrrrk 884
                     + e++pep k+kddGGslGyl+l+ll+l Glrrr+
  WP_249252303.1 852 KPEPKPEP-KEKDDGGSLGYLGLMLLGLAGLRRRN 885
                     *9999998.8899********************95 PP

>> WP_338726824.1  extracellular exonuclease ExeM [Shewanella baltica]
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ! 1609.1   7.1         0         0       1     884 []       1     869 [.       1     869 [. 0.99

  Alignments for each domain:
  == domain 1  score: 1609.1 bits;  conditional E-value: 0
  exonuc_ExeM-GG   1 mdnvkklsvlalsiaaalpllasadvliseyveGssnnkaielynsgdaavdltGyklvrykdGataasdlvaldgqtlaaksvkvivhnsaait 95 
                     mdnv+kl++++l++aaalp++asadv+i+eyveGssnnkaielynsgd+a+dl+GyklvrykdGat asd+valdgq++a k++kvi+++sa it
  WP_338726824.1   1 MDNVNKLTAISLAVAAALPMMASADVMITEYVEGSSNNKAIELYNSGDTAIDLAGYKLVRYKDGATVASDMVALDGQSIAPKTTKVILNSSAVIT 95 
                     9********************************************************************************************** PP

  exonuc_ExeM-GG  96 ldagvdsltgnlyfnGddavallkdgvvvdvvGeiptpddwakdvtlqrkadvlaastvydesewttkakdtfsGlGarsgdfvggaatepevpa 190
                     ld+gvds +g+l+fnG+daval+kd++vvd++G++ptp++w++dvtl+rk d+l+a+tv+++++w++++kdtfsGlG+++      +++epevp 
  WP_338726824.1  96 LDQGVDSYSGSLSFNGGDAVALVKDDAVVDIIGDVPTPTGWGFDVTLKRKLDALVANTVFNAAQWEQLPKDTFSGLGSLE------TPAEPEVPL 184
                     ********************************************************************************......789****** PP

  exonuc_ExeM-GG 191 fscvGakltpiydvqGaGessplvaegkyesdeevtvkGvvtargdslfkGfylqevkGdnspytsdGifvylgeaapeaiqpGvevcvqGkvke 285
                     fsc+Gak++piy+vqGaGessp v+eg++es++evtv+G+vtarg+slfkGfylqevkGdnspytsdG+fv+lgea peaiqpGvevcvqGkvke
  WP_338726824.1 185 FSCSGAKIVPIYQVQGAGESSPYVPEGAFESEAEVTVRGIVTARGESLFKGFYLQEVKGDNSPYTSDGVFVFLGEAVPEAIQPGVEVCVQGKVKE 279
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 286 yygltqidikadkklevgakaevpaaealyvadGetladaleryeGmkvkldagsdmkvtrtfsydyasrrnnlvlshkaplmkptqvyaalsde 380
                     y+gltqidikadkk+evgak+evp+a+++yvadGetla+aleryeGm+v ldagsd+k++rtfsydya+rrnn+++s+k plmk+tq+y+als e
  WP_338726824.1 280 YFGLTQIDIKADKKFEVGAKGEVPVAAPFYVADGETLAQALERYEGMNVALDAGSDLKISRTFSYDYAGRRNNMLVSYKEPLMKSTQLYPALSAE 374
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 381 aialekknkanqlfiesdlkaknGevpylpdfnaetGyirvGdqltnleGvvgysygayrlvatneitagdlirendrtdapevatkGdirvasf 475
                     a al k+n +nqlfiesd+k ++G++py+pdfn+etGyirvGdqltnl+Gv+gysygayrlvatn+itagd+ir +drtdap+vatkGd+rvasf
  WP_338726824.1 375 ATALVKSNLENQLFIESDYKPADGVIPYFPDFNVETGYIRVGDQLTNLQGVIGYSYGAYRLVATNTITAGDFIRGDDRTDAPSVATKGDLRVASF 469
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 476 nvlnlftsddavggkddadadeasrgcgrGaltleeyllqrtkivnalkemnadivGlmeiennGfgeksaiqylldalnaelpaadaykfieia 570
                     nvln+f+  d+ gg  d+++++++    rGaltlee++lqrtkiv+a+++mnadivGlmei+nnGfgeksai++l+daln++++ ++ay+f+ei+
  WP_338726824.1 470 NVLNFFN--DVDGG--DTNPSGSN----RGALTLEEMVLQRTKIVSAITAMNADIVGLMEIANNGFGEKSAIKNLVDALNEKQTPENAYSFVEIT 556
                     *******..88888..99999888....******************************************************************* PP

  exonuc_ExeM-GG 571 dadkyegkfiGsdaitvGllyraakvkpegdafvietpeqhaaeGvatreeedkvetspaydkyqrhslaqtfkikdekltvvvnhlkskGsgcl 665
                     dadky+gk++G+daitvG+lyr +kv+++g+a+ i+tpeqha++G++tr++++k+et+p++d+yqrhslaqtfki+de ltvvvnhlkskGsgcl
  WP_338726824.1 557 DADKYDGKYFGTDAITVGMLYRGGKVTLAGAAQAIDTPEQHASAGSVTRTKDGKTETNPGNDAYQRHSLAQTFKIHDESLTVVVNHLKSKGSGCL 651
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 666 edwvefeeskdpadlqGkcnefrvsaakvlGeslkdveGdllviGdlnayGledpvrvltdydaatserdivtaswttleGkvyereGskiekGy 760
                     edw++fees dpad qGkcn+frvsaakvlGe+lkdv+Gdll+iGd+nayG+edp+rvltd+da++s+rdi+taswttl+Gkv+er+GskiekGy
  WP_338726824.1 652 EDWANFEESVDPADQQGKCNAFRVSAAKVLGETLKDVKGDLLIIGDMNAYGMEDPIRVLTDFDASKSDRDIMTASWTTLDGKVFERQGSKIEKGY 746
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 761 GlinlntqahGaatysysynGelGnldhalanaslaarlvdiedwhinsvesnlfeysskytGdleksenafsasdhdpvivalsypapveevep 855
                     Glinlnt+ahGa tysysynGelGnldhalanasla+rlvdiedwhinsvesnlfey++k++Gdl+ksenafsasdhdpvivalsypapv++++p
  WP_338726824.1 747 GLINLNTKAHGAGTYSYSYNGELGNLDHALANASLAKRLVDIEDWHINSVESNLFEYGKKFSGDLAKSENAFSASDHDPVIVALSYPAPVVPPKP 841
                     *********************************************************************************************** PP

  exonuc_ExeM-GG 856 epekkkddGGslGylalallsllGlrrrk 884
                     ep++ kddGG+lGyl+lal+sl+Gl+rr+
  WP_338726824.1 842 EPTP-KDDGGALGYLGLALMSLFGLQRRR 869
                     *965.7899******************95 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                            1  (884 nodes)
Target sequences:                         12  (4815 residues searched)
Passed MSV filter:                         2  (0.166667); expected 0.2 (0.02)
Passed bias filter:                        2  (0.166667); expected 0.2 (0.02)
Passed Vit filter:                         2  (0.166667); expected 0.0 (0.001)
Passed Fwd filter:                         2  (0.166667); expected 0.0 (1e-05)
Initial search space (Z):           61295632  [as set by --Z on cmdline]
Domain search space  (domZ):               2  [number of targets reported over threshold]
# CPU time: 0.79u 0.02s 00:00:00.81 Elapsed: 00:00:00.80
# Mc/sec: 5.29
//
Query:       ExeM_NucH_DNase  [M=546]
Accession:   NF033681.1
Description: ExeM/NucH family extracellular endonuclease
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence       Description
    ------- ------ -----    ------- ------ -----   ---- --  --------       -----------
   9.1e-181  614.0   0.7   1.1e-180  613.7   0.7    1.1  1  WP_249252303.1  extracellular exonuclease ExeM [Shewanella aq
   4.5e-180  611.7   0.8   5.4e-180  611.4   0.8    1.1  1  WP_338726824.1  extracellular exonuclease ExeM [Shewanella ba


Domain annotation for each sequence (and alignments):
>> WP_249252303.1  extracellular exonuclease ExeM [Shewanella aquimarina]
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  613.7   0.7  3.6e-188  1.1e-180       1     545 [.     224     847 ..     224     848 .. 0.95

  Alignments for each domain:
  == domain 1  score: 613.7 bits;  conditional E-value: 3.6e-188
  ExeM_NucH_DNase   1 vegvVtadyqasegglkGfylqeedtdgdpatSeglFVytgsa.akevkvGdeVrvtGkvkeyygltqls..avsavevlaegeavtpvelelp 91
                      ++gvVta   + ++ +kGfylqe ++dg+p tS+g+FVy g+a  +e+++G eV+v+G vkeyygltql+  a +++ev ++  a++++++ +
   WP_249252303.1 224 IRGVVTA---RGDSLFKGFYLQEVQGDGSPYTSDGIFVYLGEAaGEEIQPGVEVCVQGLVKEYYGLTQLDikADKKMEVGENLGAPAATPFYVA 314
                      68*****...6889***************************99899************************95566677766669******8888 PP

  ExeM_NucH_DNase  92 adeear..lerlegmlvkle..qdltvtenysl...grygelvlsagerllqptqvaapgsaeaaalaaanaarrivlddgskaqnpapvpyap 178
                        e+    ler+egm v l+  +d++v++++s+   +r+++++ls++++l++ptq++++ ++ea al+++n+a++++++++ ka+n ++vpy+p
   WP_249252303.1 315 DGENLGdaLERFEGMNVVLDagSDMKVSRTFSYdyaARRNNMMLSHKAPLMKPTQLYPAMTDEAIALEKKNRANQLFVESDFKAKN-GEVPYFP 407
                      888876689***********99***********9999*************************************************.******* PP

  ExeM_NucH_DNase 179 elsadnt.lrvgdtveglegvleysfgayrlqptneaeapefvaanprtaapea.vggdlrvasfNVlNYFts............dreggetak 258
                      +++a++  +rvgd++++legv++ys+++yrl +tne+ a +f+++++r +ape  ++gd+rvasfNVlN Fts            +++ +++a+
   WP_249252303.1 408 TFNAETGyIRVGDQLTNLEGVISYSYNEYRLVATNEIVAGDFIREDDRVDAPEVaTKGDIRVASFNVLNLFTSdsdvggalnascKDQADADAS 501
                      **998777*********************************************9899****************888888877765333333799 PP

  ExeM_NucH_DNase 259 kgtnRGAksaeeferQqaKivaAlnaldAdvvgLmEienngfgedsAiadLvdaLnaaag.aetyayve.....speeeklgtDaIrvaliYkp 346
                      +g+ RGA++ e++  Q++Kiv+Al++++Ad++gLmEienngfgeds+i++Lv+ Lna++  a+ y+++e     + + e++g DaI+v+++Y+p
   WP_249252303.1 502 RGCGRGAHTLEDYLLQRTKIVSALTEMNADIIGLMEIENNGFGEDSSIQYLVNDLNANLPaADAYKFIEvadadKYKGEYVGGDAITVGMLYRP 595
                      ***********************************************************999*******99999999***************** PP

  ExeM_NucH_DNase 347 akVkpvgeaavled.....................eafdkkaReplaqtFkakesgekftvvvnHfKSKgsaca......egeadqgdgqgnwn 413
                      +kV+p+g+a v+e+                     +afdk++R++l q+F+   +gek+tvvvnH+KSKgs+c       e+++d++d qg++n
   WP_249252303.1 596 SKVEPTGDAFVIETpeqhamegvvtrgegdklelsPAFDKYQRHSLGQAFSV--NGEKLTVVVNHLKSKGSGCIedwanfEDSSDPADLQGKCN 687
                      **************************************************65..89*****************98888765556********** PP

  ExeM_NucH_DNase 414 alRvsaAkalaewlaekkkvadkdvlllGDlNaYakEdPirvLtda.................................Gytnlaeklkgeeay 474
                      ++RvsaAk +++ l+  +   ++d+l++GD+NaY+ EdPirvLtd+                                 G++nl+++++g+++y
   WP_249252303.1 688 EFRVSAAKVIGDALQ--D--VEGDLLVIGDMNAYGLEDPIRVLTDYdaansdrdimtaswttlggkvyeregskiekgyGMVNLNTQVHGTATY 777
                      ************555..5..56************************************************************************ PP

  ExeM_NucH_DNase 475 sYvfsgelGsLDhaLaseslakkvvgaeeWhiNadestaldYseeyksaedlyk.edpyrsSDHDPvivgLk 545
                      sY+++gelG+LDhaLa++sla++vv++++WhiN+ es++++Ys++y++  +l k e+++++SDHDPviv+L+
   WP_249252303.1 778 SYSYNGELGNLDHALANTSLAARVVDIQDWHINSVESNLFEYSSKYTG--ELEKsENAFSASDHDPVIVALS 847
                      ************************************************..88877***************97 PP

>> WP_338726824.1  extracellular exonuclease ExeM [Shewanella baltica]
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  611.4   0.8  1.8e-187  5.4e-180       1     545 [.     221     831 ..     221     832 .. 0.96

  Alignments for each domain:
  == domain 1  score: 611.4 bits;  conditional E-value: 1.8e-187
  ExeM_NucH_DNase   1 vegvVtadyqasegglkGfylqeedtdgdpatSeglFVytgsa.akevkvGdeVrvtGkvkeyygltqls..avsavevlaegeavtpvelelp 91
                      v+g+Vta   + e+ +kGfylqe ++d++p tS+g+FV+ g+a  +++++G eV+v+Gkvkey+gltq++  a +++ev a+ge + ++++ +
   WP_338726824.1 221 VRGIVTA---RGESLFKGFYLQEVKGDNSPYTSDGVFVFLGEAvPEAIQPGVEVCVQGKVKEYFGLTQIDikADKKFEVGAKGEVPVAAPFYVA 311
                      79*****...6789*****************************999************************977788999999999999999998 PP

  ExeM_NucH_DNase  92 adeear..lerlegmlvkle..qdltvtenysl...grygelvlsagerllqptqvaapgsaeaaalaaanaarrivlddgskaqnpapvpyap 178
                        e+ +  ler+egm v+l+  +dl++++++s+   gr++++ +s +e+l++ tq++++ saea+al ++n ++++++++++k ++ +++py+p
   WP_338726824.1 312 DGETLAqaLERYEGMNVALDagSDLKISRTFSYdyaGRRNNMLVSYKEPLMKSTQLYPALSAEATALVKSNLENQLFIESDYKPAD-GVIPYFP 404
                      888877789***********99***********999*************************************************9.******* PP

  ExeM_NucH_DNase 179 elsadnt.lrvgdtveglegvleysfgayrlqptneaeapefvaanprtaapea.vggdlrvasfNVlNYFtsdreggetakkgtnRGAksaee 270
                      +++ ++  +rvgd++++l+gv++ys+gayrl +tn+++a +f++ ++rt+ap+  ++gdlrvasfNVlN+F++   gg+t+++g+nRGA + ee
   WP_338726824.1 405 DFNVETGyIRVGDQLTNLQGVIGYSYGAYRLVATNTITAGDFIRGDDRTDAPSVaTKGDLRVASFNVLNFFND-VDGGDTNPSGSNRGALTLEE 497
                      **888767********************************************99899**************98.778888************** PP

  ExeM_NucH_DNase 271 ferQqaKivaAlnaldAdvvgLmEienngfgedsAiadLvdaLnaaag.aetyayve.....speeeklgtDaIrvaliYkpakVkpvgeaavl 358
                      +  Q++Kiv+A++a++Ad+vgLmEi nngfge+sAi++LvdaLn++++ ++ y++ve     + + +++gtDaI+v+++Y+  kV+++g+a+++
   WP_338726824.1 498 MVLQRTKIVSAITAMNADIVGLMEIANNGFGEKSAIKNLVDALNEKQTpENAYSFVEitdadKYDGKYFGTDAITVGMLYRGGKVTLAGAAQAI 591
                      *********************************************9999*********99999999**************************** PP

  ExeM_NucH_DNase 359 ed.....................eafdkkaReplaqtFkakesgekftvvvnHfKSKgsaca......egeadqgdgqgnwnalRvsaAkalae 425
                      ++                     +  d ++R++laqtFk   ++e++tvvvnH+KSKgs+c       e++ d++d+qg++na+RvsaAk l+e
   WP_338726824.1 592 DTpeqhasagsvtrtkdgktetnPGNDAYQRHSLAQTFK--IHDESLTVVVNHLKSKGSGCLedwanfEESVDPADQQGKCNAFRVSAAKVLGE 683
                      ********************998888888*********5..589*****************97777654445********************** PP

  ExeM_NucH_DNase 426 wlaekkkvadkdvlllGDlNaYakEdPirvLtda.................................GytnlaeklkgeeaysYvfsgelGsLD 486
                       l++     ++d+l++GD+NaY++EdPirvLtd                                  G++nl++k++g+ +ysY+++gelG+LD
   WP_338726824.1 684 TLKDV----KGDLLIIGDMNAYGMEDPIRVLTDFdasksdrdimtaswttldgkvferqgskiekgyGLINLNTKAHGAGTYSYSYNGELGNLD 773
                      66555....5***********************99*********************************************************** PP

  ExeM_NucH_DNase 487 haLaseslakkvvgaeeWhiNadestaldYseeyksaedlyk.edpyrsSDHDPvivgLk 545
                      haLa++slak++v++e+WhiN+ es++++Y +++++  dl+k e+++++SDHDPviv+L+
   WP_338726824.1 774 HALANASLAKRLVDIEDWHINSVESNLFEYGKKFSG--DLAKsENAFSASDHDPVIVALS 831
                      ************************************..88888***************97 PP

(...)
```

```json
{
  "WP_249252303.1": {
    "NF033680": {
      "accession": "NF033680",
      "name": "exonuc_ExeM-GG",
      "description": "extracellular exonuclease ExeM",
      "evalue": "0",
      "score": "1609.3",
      "qlen": "884",
      "bias": "7.1",
      "member_db": "ncbifam",
      "version": "14.0",
      "model-ac": "NF033680",
      "locations": [
        {
          "start": "1",
          "end": "885",
          "representative": "",
          "hmmStart": "1",
          "hmmEnd": "884",
          "hmmLength": "884",
          "rawHmmBounds": "[]",
          "hmmBounds": "Complete",
          "evalue": "0",
          "score": "1701.2",
          "envelopeStart": "1",
          "envelopeEnd": "885",
          "postProcessed": "false",
          "alignment": "MDNVKKLSVLALSIAAALPVSAKADLLISEYVEGNSNNKAIELYNSGDAALDLTGYKLVRYKDGATNALDMVALDGLTIDAKKLKVIVHPSAEIALGSDVDSSTGNLYFNGGDAVALVKDGAVVDIIGAVPTPSGWGMNTTMQRNTDALSASNVFVESQWKVLPVDTFSGLGNLDG---SSAPTEPEAPAFTCTGAKITPIYSVQGAGKSSPLIADGEYASKEEVTIRGVVTARGDSLFKGFYLQEVQGDGSPYTSDGIFVYLGEAAGEEIQPGVEVCVQGLVKEYYGLTQLDIKADKKMEVGENLGAPAATPFYVADGENLGDALERFEGMNVVLDAGSDMKVSRTFSYDYAARRNNMMLSHKAPLMKPTQLYPAMTDEAIALEKKNRANQLFVESDFKAKNGEVPYFPTFNAETGYIRVGDQLTNLEGVISYSYNEYRLVATNEIVAGDFIREDDRVDAPEVATKGDIRVASFNVLNLFTSDSDVGGalnascKDQADAD-ASRGCGRGAHTLEDYLLQRTKIVSALTEMNADIIGLMEIENNGFGEDSSIQYLVNDLNANLPAADAYKFIEVADADKYKGEYVGGDAITVGMLYRPSKVEPTGDAFVIETPEQHAMEGVVTRGEGDKLELSPAFDKYQRHSLGQAFSVNGEKLTVVVNHLKSKGSGCIEDWANFEDSSDPADLQGKCNEFRVSAAKVIGDALQDVEGDLLVIGDMNAYGLEDPIRVLTDYDAANSDRDIMTASWTTLGGKVYEREGSKIEKGYGMVNLNTQVHGTATYSYSYNGELGNLDHALANTSLAARVVDIQDWHINSVESNLFEYSSKYTGELEKSENAFSASDHDPVIVALSYPAPKPEPKPEP-KEKDDGGSLGYLGLMLLGLAGLRRRN",
          "cigar_alignment": "176M3D310M6I7M1D360M1D26M"
        },
        {
          "start": "1",
          "end": "869",
          "representative": "",
          "hmmStart": "1",
          "hmmEnd": "884",
          "hmmLength": "884",
          "rawHmmBounds": "[]",
          "hmmBounds": "Complete",
          "evalue": "0",
          "score": "1609.1",
          "envelopeStart": "1",
          "envelopeEnd": "869",
          "postProcessed": "false",
          "alignment": "MDNVNKLTAISLAVAAALPMMASADVMITEYVEGSSNNKAIELYNSGDTAIDLAGYKLVRYKDGATVASDMVALDGQSIAPKTTKVILNSSAVITLDQGVDSYSGSLSFNGGDAVALVKDDAVVDIIGDVPTPTGWGFDVTLKRKLDALVANTVFNAAQWEQLPKDTFSGLGSLE------TPAEPEVPLFSCSGAKIVPIYQVQGAGESSPYVPEGAFESEAEVTVRGIVTARGESLFKGFYLQEVKGDNSPYTSDGVFVFLGEAVPEAIQPGVEVCVQGKVKEYFGLTQIDIKADKKFEVGAKGEVPVAAPFYVADGETLAQALERYEGMNVALDAGSDLKISRTFSYDYAGRRNNMLVSYKEPLMKSTQLYPALSAEATALVKSNLENQLFIESDYKPADGVIPYFPDFNVETGYIRVGDQLTNLQGVIGYSYGAYRLVATNTITAGDFIRGDDRTDAPSVATKGDLRVASFNVLNFFN--DVDGG--DTNPSGSN----RGALTLEEMVLQRTKIVSAITAMNADIVGLMEIANNGFGEKSAIKNLVDALNEKQTPENAYSFVEITDADKYDGKYFGTDAITVGMLYRGGKVTLAGAAQAIDTPEQHASAGSVTRTKDGKTETNPGNDAYQRHSLAQTFKIHDESLTVVVNHLKSKGSGCLEDWANFEESVDPADQQGKCNAFRVSAAKVLGETLKDVKGDLLIIGDMNAYGMEDPIRVLTDFDASKSDRDIMTASWTTLDGKVFERQGSKIEKGYGLINLNTKAHGAGTYSYSYNGELGNLDHALANASLAKRLVDIEDWHINSVESNLFEYGKKFSGDLAKSENAFSASDHDPVIVALSYPAPVVPPKPEPTP-KDDGGALGYLGLALMSLFGLQRRR",
          "cigar_alignment": "175M6D301M2D5M2D8M4D356M1D24M"
        }
      ]
    },
    "NF033681": {
      "accession": "NF033681",
      "name": "ExeM_NucH_DNase",
      "description": "ExeM/NucH family extracellular endonuclease",
      "evalue": "4.5e-180",
      "score": "611.7",
      "qlen": "546",
      "bias": "0.8",
      "member_db": "ncbifam",
      "version": "14.0",
      "model-ac": "NF033681",
      "locations": [
        {
          "start": "224",
          "end": "847",
          "representative": "",
          "hmmStart": "1",
          "hmmEnd": "545",
          "hmmLength": "546",
          "rawHmmBounds": "[.",
          "hmmBounds": "N-terminal complete",
          "evalue": "1.1e-180",
          "score": "613.7",
          "envelopeStart": "224",
          "envelopeEnd": "848",
          "postProcessed": "false",
          "alignment": "IRGVVTA---RGDSLFKGFYLQEVQGDGSPYTSDGIFVYLGEAaGEEIQPGVEVCVQGLVKEYYGLTQLDikADKKMEVGENLGAPAATPFYVADGENLGdaLERFEGMNVVLDagSDMKVSRTFSYdyaARRNNMMLSHKAPLMKPTQLYPAMTDEAIALEKKNRANQLFVESDFKAKN-GEVPYFPTFNAETGyIRVGDQLTNLEGVISYSYNEYRLVATNEIVAGDFIREDDRVDAPEVaTKGDIRVASFNVLNLFTSdsdvggalnascKDQADADASRGCGRGAHTLEDYLLQRTKIVSALTEMNADIIGLMEIENNGFGEDSSIQYLVNDLNANLPaADAYKFIEvadadKYKGEYVGGDAITVGMLYRPSKVEPTGDAFVIETpeqhamegvvtrgegdklelsPAFDKYQRHSLGQAFSV--NGEKLTVVVNHLKSKGSGCIedwanfEDSSDPADLQGKCNEFRVSAAKVIGDALQ--D--VEGDLLVIGDMNAYGLEDPIRVLTDYdaansdrdimtaswttlggkvyeregskiekgyGMVNLNTQVHGTATYSYSYNGELGNLDHALANTSLAARVVDIQDWHINSVESNLFEYSSKYTG--ELEKsENAFSASDHDPVIVALS",
          "cigar_alignment": "7M3D33M1I26M2I28M2I12M2I11M3I50M1D14M1I46M1I18M12I69M1I8M5I34M21I17M2D20M6I29M2D1M2D26M33I63M2D4M1I17M"
        },
        {
          "start": "221",
          "end": "831",
          "representative": "",
          "hmmStart": "1",
          "hmmEnd": "545",
          "hmmLength": "546",
          "rawHmmBounds": "[.",
          "hmmBounds": "N-terminal complete",
          "evalue": "5.4e-180",
          "score": "611.4",
          "envelopeStart": "221",
          "envelopeEnd": "832",
          "postProcessed": "false",
          "alignment": "VRGIVTA---RGESLFKGFYLQEVKGDNSPYTSDGVFVFLGEAvPEAIQPGVEVCVQGKVKEYFGLTQIDikADKKFEVGAKGEVPVAAPFYVADGETLAqaLERYEGMNVALDagSDLKISRTFSYdyaGRRNNMLVSYKEPLMKSTQLYPALSAEATALVKSNLENQLFIESDYKPAD-GVIPYFPDFNVETGyIRVGDQLTNLQGVIGYSYGAYRLVATNTITAGDFIRGDDRTDAPSVaTKGDLRVASFNVLNFFND-VDGGDTNPSGSNRGALTLEEMVLQRTKIVSAITAMNADIVGLMEIANNGFGEKSAIKNLVDALNEKQTpENAYSFVEitdadKYDGKYFGTDAITVGMLYRGGKVTLAGAAQAIDTpeqhasagsvtrtkdgktetnPGNDAYQRHSLAQTFK--IHDESLTVVVNHLKSKGSGCLedwanfEESVDPADQQGKCNAFRVSAAKVLGETLKDV----KGDLLIIGDMNAYGMEDPIRVLTDFdasksdrdimtaswttldgkvferqgskiekgyGLINLNTKAHGAGTYSYSYNGELGNLDHALANASLAKRLVDIEDWHINSVESNLFEYGKKFSG--DLAKsENAFSASDHDPVIVALS",
          "cigar_alignment": "7M3D33M1I26M2I28M2I12M2I11M3I50M1D14M1I46M1I18M1D68M1I8M5I34M21I16M2D21M6I31M4D25M33I63M2D4M1I17M"
        }
      ]
    }
  },
  "WP_338726824.1": {
    "NF033680": {
      "accession": "NF033680",
      "name": "exonuc_ExeM-GG",
      "description": "extracellular exonuclease ExeM",
      "evalue": "0",
      "score": "1609.3",
      "qlen": "884",
      "bias": "7.1",
      "member_db": "ncbifam",
      "version": "14.0",
      "model-ac": "NF033680",
      "locations": [
        {
          "start": "1",
          "end": "885",
          "representative": "",
          "hmmStart": "1",
          "hmmEnd": "884",
          "hmmLength": "884",
          "rawHmmBounds": "[]",
          "hmmBounds": "Complete",
          "evalue": "0",
          "score": "1701.2",
          "envelopeStart": "1",
          "envelopeEnd": "885",
          "postProcessed": "false",
          "alignment": "MDNVKKLSVLALSIAAALPVSAKADLLISEYVEGNSNNKAIELYNSGDAALDLTGYKLVRYKDGATNALDMVALDGLTIDAKKLKVIVHPSAEIALGSDVDSSTGNLYFNGGDAVALVKDGAVVDIIGAVPTPSGWGMNTTMQRNTDALSASNVFVESQWKVLPVDTFSGLGNLDG---SSAPTEPEAPAFTCTGAKITPIYSVQGAGKSSPLIADGEYASKEEVTIRGVVTARGDSLFKGFYLQEVQGDGSPYTSDGIFVYLGEAAGEEIQPGVEVCVQGLVKEYYGLTQLDIKADKKMEVGENLGAPAATPFYVADGENLGDALERFEGMNVVLDAGSDMKVSRTFSYDYAARRNNMMLSHKAPLMKPTQLYPAMTDEAIALEKKNRANQLFVESDFKAKNGEVPYFPTFNAETGYIRVGDQLTNLEGVISYSYNEYRLVATNEIVAGDFIREDDRVDAPEVATKGDIRVASFNVLNLFTSDSDVGGalnascKDQADAD-ASRGCGRGAHTLEDYLLQRTKIVSALTEMNADIIGLMEIENNGFGEDSSIQYLVNDLNANLPAADAYKFIEVADADKYKGEYVGGDAITVGMLYRPSKVEPTGDAFVIETPEQHAMEGVVTRGEGDKLELSPAFDKYQRHSLGQAFSVNGEKLTVVVNHLKSKGSGCIEDWANFEDSSDPADLQGKCNEFRVSAAKVIGDALQDVEGDLLVIGDMNAYGLEDPIRVLTDYDAANSDRDIMTASWTTLGGKVYEREGSKIEKGYGMVNLNTQVHGTATYSYSYNGELGNLDHALANTSLAARVVDIQDWHINSVESNLFEYSSKYTGELEKSENAFSASDHDPVIVALSYPAPKPEPKPEP-KEKDDGGSLGYLGLMLLGLAGLRRRN",
          "cigar_alignment": "176M3D310M6I7M1D360M1D26M"
        },
        {
          "start": "1",
          "end": "869",
          "representative": "",
          "hmmStart": "1",
          "hmmEnd": "884",
          "hmmLength": "884",
          "rawHmmBounds": "[]",
          "hmmBounds": "Complete",
          "evalue": "0",
          "score": "1609.1",
          "envelopeStart": "1",
          "envelopeEnd": "869",
          "postProcessed": "false",
          "alignment": "MDNVNKLTAISLAVAAALPMMASADVMITEYVEGSSNNKAIELYNSGDTAIDLAGYKLVRYKDGATVASDMVALDGQSIAPKTTKVILNSSAVITLDQGVDSYSGSLSFNGGDAVALVKDDAVVDIIGDVPTPTGWGFDVTLKRKLDALVANTVFNAAQWEQLPKDTFSGLGSLE------TPAEPEVPLFSCSGAKIVPIYQVQGAGESSPYVPEGAFESEAEVTVRGIVTARGESLFKGFYLQEVKGDNSPYTSDGVFVFLGEAVPEAIQPGVEVCVQGKVKEYFGLTQIDIKADKKFEVGAKGEVPVAAPFYVADGETLAQALERYEGMNVALDAGSDLKISRTFSYDYAGRRNNMLVSYKEPLMKSTQLYPALSAEATALVKSNLENQLFIESDYKPADGVIPYFPDFNVETGYIRVGDQLTNLQGVIGYSYGAYRLVATNTITAGDFIRGDDRTDAPSVATKGDLRVASFNVLNFFN--DVDGG--DTNPSGSN----RGALTLEEMVLQRTKIVSAITAMNADIVGLMEIANNGFGEKSAIKNLVDALNEKQTPENAYSFVEITDADKYDGKYFGTDAITVGMLYRGGKVTLAGAAQAIDTPEQHASAGSVTRTKDGKTETNPGNDAYQRHSLAQTFKIHDESLTVVVNHLKSKGSGCLEDWANFEESVDPADQQGKCNAFRVSAAKVLGETLKDVKGDLLIIGDMNAYGMEDPIRVLTDFDASKSDRDIMTASWTTLDGKVFERQGSKIEKGYGLINLNTKAHGAGTYSYSYNGELGNLDHALANASLAKRLVDIEDWHINSVESNLFEYGKKFSGDLAKSENAFSASDHDPVIVALSYPAPVVPPKPEPTP-KDDGGALGYLGLALMSLFGLQRRR",
          "cigar_alignment": "175M6D301M2D5M2D8M4D356M1D24M"
        }
      ]
    },
    "NF033681": {
      "accession": "NF033681",
      "name": "ExeM_NucH_DNase",
      "description": "ExeM/NucH family extracellular endonuclease",
      "evalue": "4.5e-180",
      "score": "611.7",
      "qlen": "546",
      "bias": "0.8",
      "member_db": "ncbifam",
      "version": "14.0",
      "model-ac": "NF033681",
      "locations": [
        {
          "start": "224",
          "end": "847",
          "representative": "",
          "hmmStart": "1",
          "hmmEnd": "545",
          "hmmLength": "546",
          "rawHmmBounds": "[.",
          "hmmBounds": "N-terminal complete",
          "evalue": "1.1e-180",
          "score": "613.7",
          "envelopeStart": "224",
          "envelopeEnd": "848",
          "postProcessed": "false",
          "alignment": "IRGVVTA---RGDSLFKGFYLQEVQGDGSPYTSDGIFVYLGEAaGEEIQPGVEVCVQGLVKEYYGLTQLDikADKKMEVGENLGAPAATPFYVADGENLGdaLERFEGMNVVLDagSDMKVSRTFSYdyaARRNNMMLSHKAPLMKPTQLYPAMTDEAIALEKKNRANQLFVESDFKAKN-GEVPYFPTFNAETGyIRVGDQLTNLEGVISYSYNEYRLVATNEIVAGDFIREDDRVDAPEVaTKGDIRVASFNVLNLFTSdsdvggalnascKDQADADASRGCGRGAHTLEDYLLQRTKIVSALTEMNADIIGLMEIENNGFGEDSSIQYLVNDLNANLPaADAYKFIEvadadKYKGEYVGGDAITVGMLYRPSKVEPTGDAFVIETpeqhamegvvtrgegdklelsPAFDKYQRHSLGQAFSV--NGEKLTVVVNHLKSKGSGCIedwanfEDSSDPADLQGKCNEFRVSAAKVIGDALQ--D--VEGDLLVIGDMNAYGLEDPIRVLTDYdaansdrdimtaswttlggkvyeregskiekgyGMVNLNTQVHGTATYSYSYNGELGNLDHALANTSLAARVVDIQDWHINSVESNLFEYSSKYTG--ELEKsENAFSASDHDPVIVALS",
          "cigar_alignment": "7M3D33M1I26M2I28M2I12M2I11M3I50M1D14M1I46M1I18M12I69M1I8M5I34M21I17M2D20M6I29M2D1M2D26M33I63M2D4M1I17M"
        },
        {
          "start": "221",
          "end": "831",
          "representative": "",
          "hmmStart": "1",
          "hmmEnd": "545",
          "hmmLength": "546",
          "rawHmmBounds": "[.",
          "hmmBounds": "N-terminal complete",
          "evalue": "5.4e-180",
          "score": "611.4",
          "envelopeStart": "221",
          "envelopeEnd": "832",
          "postProcessed": "false",
          "alignment": "VRGIVTA---RGESLFKGFYLQEVKGDNSPYTSDGVFVFLGEAvPEAIQPGVEVCVQGKVKEYFGLTQIDikADKKFEVGAKGEVPVAAPFYVADGETLAqaLERYEGMNVALDagSDLKISRTFSYdyaGRRNNMLVSYKEPLMKSTQLYPALSAEATALVKSNLENQLFIESDYKPAD-GVIPYFPDFNVETGyIRVGDQLTNLQGVIGYSYGAYRLVATNTITAGDFIRGDDRTDAPSVaTKGDLRVASFNVLNFFND-VDGGDTNPSGSNRGALTLEEMVLQRTKIVSAITAMNADIVGLMEIANNGFGEKSAIKNLVDALNEKQTpENAYSFVEitdadKYDGKYFGTDAITVGMLYRGGKVTLAGAAQAIDTpeqhasagsvtrtkdgktetnPGNDAYQRHSLAQTFK--IHDESLTVVVNHLKSKGSGCLedwanfEESVDPADQQGKCNAFRVSAAKVLGETLKDV----KGDLLIIGDMNAYGMEDPIRVLTDFdasksdrdimtaswttldgkvferqgskiekgyGLINLNTKAHGAGTYSYSYNGELGNLDHALANASLAKRLVDIEDWHINSVESNLFEYGKKFSG--DLAKsENAFSASDHDPVIVALS",
          "cigar_alignment": "7M3D33M1I26M2I28M2I12M2I11M3I50M1D14M1I46M1I18M1D68M1I8M5I34M21I16M2D21M6I31M4D25M33I63M2D4M1I17M"
        }
      ]
    },
    "TIGR03501": {
      "accession": "TIGR03501",
      "name": "GlyGly_CTERM",
      "description": "GlyGly-CTERM protein-sorting domain",
      "evalue": "20",
      "score": "17.2",
      "qlen": "22",
      "bias": "3.0",
      "member_db": "ncbifam",
      "version": "14.0",
      "model-ac": "TIGR03501",
      "locations": [
        {
          "start": "849",
          "end": "869",
          "representative": "",
          "hmmStart": "2",
          "hmmEnd": "22",
          "hmmLength": "22",
          "rawHmmBounds": ".]",
          "hmmBounds": "C-terminal complete",
          "evalue": "46",
          "score": "16.1",
          "envelopeStart": "848",
          "envelopeEnd": "869",
          "postProcessed": "false",
          "alignment": "GGALGYLGLALMSLFGLQRRR",
          "cigar_alignment": "21M"
        }
      ]
    }
  },
  "OUTROSEMLOOKUP": {
    "TIGR01274": {
      "accession": "TIGR01274",
      "name": "ACC_deam",
      "description": "1-aminocyclopropane-1-carboxylate deaminase",
      "evalue": "7.3e-197",
      "score": "664.8",
      "qlen": "337",
      "bias": "0.1",
      "member_db": "ncbifam",
      "version": "14.0",
      "model-ac": "TIGR01274",
      "locations": [
        {
          "start": "2",
          "end": "317",
          "representative": "",
          "hmmStart": "1",
          "hmmEnd": "316",
          "hmmLength": "337",
          "rawHmmBounds": "[.",
          "hmmBounds": "N-terminal complete",
          "evalue": "1.8e-185",
          "score": "627.4",
          "envelopeStart": "2",
          "envelopeEnd": "317",
          "postProcessed": "false",
          "alignment": "NLNRFERYPLTFGPSPITPLKRLSQHLGGKVELYAKREDCNSGLAFGGNKTRKLEYLIPEAIEQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKCVLVQENWVNYSDAVYDRVGNIEMSRIMGADVRLDAAGFDIGIRPSWEKAMSDVVEQGGKPFPIPAGCSEHPYGGLGFVGFAEEVRQQEKELGFKFDYIVVCSVTGSTQAGMVVGFAADGRSKNVIGIDASAKPEQTKAQILRIARHTAELVELGREITEEDVVLDTRFAYPEYGLPNEGTLEAIRLCGSLEGVLTDPVYEGKSMHGMIEMVRRGEFPEGSMI",
          "cigar_alignment": "316M"
        }
      ]
    }
  }
}
```
