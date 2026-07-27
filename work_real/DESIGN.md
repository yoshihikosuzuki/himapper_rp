# HG002 chr6 mapping benchmark

## 目的

HG002 HiFi readsをT2T-HG002 v1.1のchr6 maternal/paternal panelへmappingし、
Himapperのseed構成3通りとminimap2/winnowmap2を同一入力・同一評価条件で比較する。

現時点の本命候補は `ref-seed_read-seed` だが、比較完了前にproductionとは呼ばない。
過去の `work_real/out/` の結果は条件が揃っていないため、新しいfrom-scratch比較の
canonical resultにはしない。

## 共通mapping input

- FASTA: `reads/chr6.reads.fa`
- ordered names: `reads/chr6.reads.names`

全Himapper config、minimap2、winnowmap2が同じread集合と順序を使う。

`reads/reads.map.paf` から、次を満たすreadを選ぶ。

- targetが `chr6_MATERNAL` または `chr6_PATERNAL`
- read length >= 10 kb
- aligned query coverage >= 0.90

`02_prepare_reads.sh` は最長alignmentを選んだ後にchr6へ限定し、
`reads/chr6.reads.mapping_truth.tsv` とmapping inputを生成する。coverageを含む全mapperは
この広いread集合を使う。

`reads/chr6.reads.confident.mapping_truth.tsv` は上記mapping truthのうち
MAPQ >= 60の部分集合であり、stratificationとaccuracy評価の固定母集団に使う。

## Trio classification

`03_prepare_trio.sh` はreadをtrio binningで分類する。

- `trio/chr6.reads.hap.tsv`: trioでpaternalまたはmaternalに分類されたread

`reads/chr6.reads.confident.mapping_truth.tsv` をself-map difficulty判定と
mapping評価に使う。

## 浅いstage構成

```text
work_real/
  ref/
  reads/
  trio/
  regions/
  himapper/
    ref-all_read-seed/
    ref-seed_read-all/
    ref-seed_read-seed/
  mm2/
  wm2/
  eval/
```

`00`〜`08`が主比較pipeline、`09`は任意の追加解析である。
Himapperの各`07`は1 kb window、54–68 Mbでminimap2/winnowmap2と
比較するcoverage plotも生成する。
mapper別scriptをさらに下位directoryへ分けず、top level一段に保つ。

## Himapper config IDs

| config ID | reference seed spec | read seed spec | tier source |
|---|---|---|---|
| `ref-all_read-seed` | `all` | `.reads.minimax.kpos` | read (`CN_TIER_SRC=1`) |
| `ref-seed_read-all` | `.ref.minimax.kpos` | `all` | reference (`CN_TIER_SRC=0`) |
| `ref-seed_read-seed` | `.ref.minimax.kpos` | `.reads.minimax.kpos` | read (`CN_TIER_SRC=1`) |

config IDは、比較で変えるreference/read両側のseed specを機械的に表す安定IDである。
「採用済み」「本番」等の状態は含めず、script名、output directory、log、評価labelを対応させる。

`ref-seed_read-all` では `w0` を、copy-aware経路を有効にしたままbonusを0にする
no-copy-bonus baselineとして扱う。W=11.5候補とは別設定である。

`minimaxmer` はwhole DBとwhole FastK own/cross profilesを入力し、`-T`でDB blockを
並列処理する。referenceは `.ref.minimax.kpos`、readsは
`.reads.minimax.kpos` という単一containerと、それぞれのblock境界を保持する
`.kpos.idx` sidecarを生成する。entryのIDはblock-localで、各block内はcode順に
sortされる。3構成はすべて`DBsplit -s200`したwhole reference/read DBを
`himapper`へ1回だけ渡し、reads blockの反復とreference KPOSの再利用は同一process内で
行う。

## Dataとoutputの配置

```text
work_real/
  ref/
  reads/
  trio/
  regions/
  himapper/
    ref-all_read-seed/
    ref-seed_read-all/
    ref-seed_read-seed/
  mm2/
  wm2/
  eval/
```

- DAM、DBsplit metadata、FastK own/cross profilesはHimapper source変更に対する共通入力なので
  `himapper/` に置く。ただしmapping reads、reference、K、FastK条件には依存する。
- minimax `.kpos` はminimaxmerの実装・parameterに依存するため
  `himapper/` に置く。reference/readそれぞれをhidden containerとsidecar
  `.ref.minimax.kpos{,.idx}` / `.reads.minimax.kpos{,.idx}` として置く。
- mapper結果directoryの中へ別configが依存するshared inputを置かない。
- 旧 `himapper/ref-all_read-seed/work/` はlegacy artifactであり、新pipelineからは参照しない。

mapping readsの集合と順序は、生成scriptとmanifest/checksumで検証する。

## Stratification derived products

`01_prepare_reference.sh` が `ref/variants.tsv` を生成する。

`03_prepare_trio.sh` が `trio/chr6.hapmer.tsv` を生成する。

`04_prepare_eval.sh` が次を生成する。

- `eval/alignment_strata.tsv`: mapping-truth spanのminimap2 self-map difficulty
  proxy。`no_map`、`dup`、`hard`、`med`、`easy`、NM=0 alignment数、
  最も近い別locusのNM
- `eval/hapmer_strata.tsv`: hapmer countとtrio classificationによるstratum
- `eval/variant_strata.tsv`: locus内のMAT/PAT variantの有無とvariant数

accuracy評価の母集団はMAPQ >= 60のmapping truth全体に固定する。
alignment strataやvariantの有無によるread除外は行わず、`eval.py --strata`
の`alignment`、`hapmer`、`variant`から層別軸を選ぶ。`08_evaluate.sh`は
overallを一度だけ計算し、hapmerとalignmentの両stratificationを同じreportへ
出力する。

downloadしたBED等のsource assetそのものはここでいう派生物に含めない。

## Mapping outputs

mapping outputはmapperを実行した結果・log・timingであり、mapping truth、
DAM/FastK、seed cacheは含めない。

- Himapper (`ref-all_read-seed`, `ref-seed_read-all`, `ref-seed_read-seed`): `himapper/<config ID>/mapping.tsv`, `map*.time`, `coverage.{png,pdf}`, retained block LAS (`reads.N.ref*.las`)
- minimap2: `mm2/chr6.mmi`, `index.time`, `map.paf`, `mapping.tsv`, `map.time`, sorted BAM, `bam.time`
- winnowmap2: `wm2/chr6.meryl`, repetitive k-mers, preparation times, `map.paf`, `mapping.tsv`, `map.time`, sorted BAM, `bam.time`

各commandのstdout/stderrは独立logへ分けず、`#SBATCH -o`のstage logへ集約する。

minimap2/winnowmap2の `map.time` は `-c` によるbase-level alignmentとPAF出力を測る。
`-a | samtools sort` はmapping+sort比較として別実行し、`bam.time` に記録する。
BAM index作成は含めない。

正規化した`mapping.tsv`をprimary accuracy評価とcoverage plotの両方に使う。
`eval/stratified.txt` は過去分を削除済みで、`08_evaluate.sh` を実行した新runだけが再生成する。

## 実行順

全scriptは `work_real/` をcurrent directoryとして投入または実行する。

```bash
sbatch 00_download.sh
sbatch 01_prepare_reference.sh
sbatch 02_prepare_reads.sh
sbatch 03_prepare_trio.sh
sbatch 04_prepare_eval.sh
# 04完了後
sbatch 05_prepare_himapper.sh
sbatch 06_map_minimap2.sh
sbatch 06_map_winnowmap2.sh

# DAM/profileを再利用してKPOSだけ更新する場合
sbatch --export=ALL,SEEDS_ONLY=1 05_prepare_himapper.sh

# 05と両方の06完了後
sbatch 07_map_himapper_ref-seed_read-seed.sh
# 次の2構成は追加比較を行う場合に実行
sbatch 07_map_himapper_ref-all_read-seed.sh
sbatch 07_map_himapper_ref-seed_read-all.sh

# 必要なmapping構成の完了後
sbatch 08_evaluate.sh
```

この文書は設計と実行順を示す。
