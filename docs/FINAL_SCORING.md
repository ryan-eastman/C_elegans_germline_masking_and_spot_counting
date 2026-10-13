# Final scoring: the one-time test scoring of the RAD-51 detectors

Written 2026-10-10 (package E1). **No step of this plan has been run.** It is the plan for the one test
scoring that `docs/BENCHMARK_PROTOCOL.md` allows per method family (sections 7, 8 and 10): the decisions to
take first, the frozen methods, the exact commands in order, how the win conditions are read, and the
changelog text to append before and after. No test gonad was scored, run or looked at while it was written.

The test gonads (protocol section 3): n2rad51_HS_herm_002, n2rad51_HS_male_001, n2rad51_noHS_herm_002,
n2rad51_noHS_male_002.

## 1. Decisions the project lead takes before anything runs

### Decision 1. Wait for the re-exported heat-shock male spots, or score now

* The two heat-shock males of train, n2rad51_HS_male_000 and n2rad51_HS_male_002, were exported from Imaris
  with region growing on; their spot truth is unusable (protocol section 4.1) and no heat-shock male has
  informed any choice: not the D3 searches, not net3d's training or selection, not the leave-one-gonad-out
  estimates. Val holds hermaphrodites only.
* The deployed recipe's heat-shock undercount, the in-house detector's main target, is not visible in train
  or val: on the two heat-shock hermaphrodites the deployed SpotMAX recipe over-counts (docs/BENCHMARK.md,
  "What the heat-shock gonads show"). The only heat-shock male with clean spot truth is n2rad51_HS_male_001,
  a test gonad. Scoring now therefore tests the main claim on one gonad that no method has been shaped on.
* Waiting means: re-export both projects without region growing, record their new size and sha256 in the
  split record with a changelog line (section 4.1), rebuild their truth files, and then decide whether they
  enter training. If they do, every learned or tuned method changes (net3d retrained on eight gonads, the
  classical searches rerun or at least rescored, the leave-one-gonad-out estimates redone on eight gonads):
  a declared revision before any test scoring, which costs a day of compute and the Imaris work. If they
  only enter the leave-one-gonad-out estimate (as evidence, not as training data), nothing frozen changes.
* Scoring now means the write-up states that no heat-shock male informed any choice and that the
  heat-shock-male result rests on one test gonad.

Status (2026-10-11, protocol changelog entry 6): the project lead replied "Yes to everything and merged" to the
message that listed the five decisions. That does not settle this choice by itself. It is read as waiting
for the re-export, the reading that keeps both options open: no test gonad is scored until he confirms one
option explicitly.

Status (2026-10-12, protocol changelog entry 9): still open. The project lead accepted the 2026-01-05
heat-shocked males as training data ("2 whatever you recommend") and allowed a check of the two training males'
2025-12-31 exports ("1 yes"). Entry 9 adds three or four heat-shocked males to train (four if the xlsx spots of
20260105_n2_hs_1_MALE_1 verify), so the main reason to wait is met without the Imaris re-export. The methods this
decision names change in a revision declared there before any test scoring: net3d retrained by the frozen recipe,
the tuned log's D3 procedure replayed, leave-one-gonad-out redone. SpotMAX re-tuned and SpotMAX deployed do not
change. The check of the 2025-12-31 exports, made by name, size and file time only with no file of the test male
and no 20251105_* file opened, reported (entry 9 (i) gives it in full): which file is which is known with high
confidence, confirmed by a second agent (MALE_1 is the sealed test male n2rad51_HS_male_001, its export mirroring
the registered 10,411 spots and 68 surfaces; MALE_2 is n2rad51_HS_male_002; MALE_3 is n2rad51_HS_male_000). The
exports are separable but not clean: MALE_2 looks region-grown and holds 13,690 to 13,695 spots with every
statistic (258 fewer than the registered 13,951) or exactly 13,951 with one statistic missing, so it gives no clean
spot set with region growing off; MALE_3's .xls holds about 0.75 to 0.77 million spots, and its .xlsx, of
unexplained origin, about 7,300 to 7,700 spots with about 64 surfaces, the only 2025-12-31 file that might give
n2rad51_HS_male_000 a usable spot set, which only opening it could settle; it was not opened (the check closed
without a content check). The lab's aggregate tables that very likely compile the test male's per-nucleus rows are
sealed as well. So the revision runs with the additions only, and n2rad51_HS_male_000 and _002 stay out of every
spot metric. The options are now: (1) score after the entry-9 revision (recommended: the additions supply the
heat-shocked male training signal), or (2) wait for an Imaris re-export. No test gonad is scored until he confirms
one.

Status (2026-10-12, after the CPU truth steps of entry 9; protocol changelog entry 10): under the section 4.4
criterion entry 9 declared, one addition keeps spot truth (n2rad51_20260105a_HS_male_002); three if the project
lead accepts the post hoc reading for n2rad51_20260105a_HS_male_001 and _003 (like for like, _001 lies in the
registered range and _003 does not); 20260105_n2_hs_1_MALE_1 is nucleus truth only. "Accepted" above is entry 9
(a)'s reading of "2 whatever you recommend", not his words. The export check opened no 20251105 file; for
calibration it opened about 180 non-benchmark exports and processed tables and one lab R script. Under option (1)
the heat-shocked male training signal is therefore one male unless he keeps more. Put to him with this decision,
each answer to be recorded verbatim in a changelog entry: (a) for each of n2rad51_20260105a_HS_male_001 and _003,
keep its spot truth under the post hoc reading of entry 10 (c), or not; (b) whether 20260105_n2_hs_1_MALE_1 may
enter if a later check verifies it; (c) that 2026-01-05 then cannot serve as the session of a future test set
(protocol section 10). The steps of the entry-9 revision that use spot truth wait for (a).

Status (2026-10-12, protocol changelog entry 11): the project lead answered "9 yes" to using n2rad51_HS_male_000's
recovered spots in training, "10 sure" to keeping n2rad51_20260105a_HS_male_001's spot truth under the post hoc
like-for-like reading (_003 stays out), and "okay i will re export" to re-exporting n2rad51_HS_male_002 from Imaris
(region growing off, Quality 6.847, the A and B classes kept). n2rad51_HS_male_000's spot truth now comes from its
2025-12-31 export: 7,628 curated spots at constant diameter (region growing off), a position-exact subset of the
772,485 unthresholded candidates saved in its .ims, recorded beside the split in
`tests/golden/rad51_spot_sources.json` (an earlier curated export, not the re-export of protocol section 4.1; the
split and its roles do not change), every check of entry 11 (c) passed, section 4.4 at 0.9967. The heat-shocked
male training signal is therefore three males (n2rad51_HS_male_000, n2rad51_20260105a_HS_male_001 and _002), and
leave-one-gonad-out has nine folds, ten if the re-export passes its checks. The steps of the entry-9 revision that
use spot truth wait for that re-export, which enters by the path of protocol section 4.1 and entry 10 (g)(8), or run
without it if he says so. In the template of section 4.2 the net3d brackets then read nine or ten train and val
gonads, seven or eight of them registered and two additions. Still open: the option of this decision (no test gonad
is scored until he confirms one), and questions (b) and (c) above.

Status (2026-10-12, protocol changelog entry 12, which corrects entry 11 and the status line above):
n2rad51_HS_male_000's spot truth is post hoc under rule 3 (entry 11), like n2rad51_20260105a_HS_male_001's. The
export's 7,628 spots are at a constant estimated diameter, that is, without region-grown sizes; "region growing off"
and "curated" above are inferences (the settings the export was made with and its Quality threshold are recorded
nowhere, and the .ims was saved again on 2026-01-21, after the export). Under the criteria entry 9 (e) declared for
this file, the export fails V3 on volume and n2rad51_HS_male_000 would be nucleus truth only, with no per-spot
classes; entry 11 replaced those criteria after the file had been examined. That is put to the project lead with
this decision; until he answers, n2rad51_HS_male_000 keeps its spot truth, and if he does not confirm it, a later
entry makes it nucleus truth only and the steps of the revision that used it run again before any test scoring. The
gate of the revision replaces "or run without it if he says so" above with what was put to him ("If you can do it
today, the retraining runs once with all the heat-shocked males. If not, I go ahead without it."): the steps that
use spot truth wait for the re-export of n2rad51_HS_male_002 only if it arrives on 2026-10-12; otherwise they run
with nine folds without it, and the re-export enters afterwards as a further declared development revision before
test scoring (entry 10 (g)(8)).

Status (2026-10-12, protocol changelog entry 13): the project lead replied "Whatever you think I don't have time to
export right now." (apostrophe in ASCII) to q11 (keep n2rad51_HS_male_000 in training on the adapted checks; the
recommendation was yes) and to the timing of the re-export. Read as accepting q11, so n2rad51_HS_male_000 keeps its
spot truth under the criteria of entry 11 (post hoc, entry 12 (a)), and as saying that the re-export of
n2rad51_HS_male_002 does not arrive on 2026-10-12. The steps of the entry-9 revision that use spot truth therefore run
with nine folds (the six registered, n2rad51_HS_male_000, n2rad51_20260105a_HS_male_001 and _002); the re-export,
when it comes, enters by the path of protocol section 4.1 (entries 10 (g)(8) and 12 (d)) as a further declared
development revision before test scoring, with ten folds. In the template of section 4.2 the net3d brackets read
nine train and val gonads, seven of them registered, and two additions, until that revision. The message stated,
without asking, that 2026-01-05 cannot be the session of a future test set and that 20260105_n2_hs_1_MALE_1 could join
if its Spots are saved again in Imaris; questions (b) and (c) above stay open. The option of this decision is still
his: no test gonad is scored until he confirms one.

Status (2026-10-12, protocol changelog entry 14, which corrects entry 13 before the revision runs): "right now" is
narrower than "on 2026-10-12". The reply is read, with "Whatever you think", as accepting the course the message named
without the re-export ("Otherwise it runs tonight on the 9 I have"): the revision runs tonight on nine folds, and a
re-export that arrives later on 2026-10-12 enters as the further declared revision above, not by holding the run,
which departs from the letter of entry 12 (c). Nothing else above changes: the option of this decision and questions
(b) and (c) stay open.

Status (2026-10-13, protocol changelog entry 15): the entry-9 revision ran on nine folds on 2026-10-12, 18:07 to 23:52
MDT, from commit 08f6285, and every check passed; no test gonad was run, scored or opened. Out of sample on the three
heat-shocked males (250 truth nuclei, about 17 truth foci per traced nucleus), truth-mask CCC and bias in foci per truth
nucleus: net3d 0.621 and +3.42, the tuned log 0.635 and -1.84, SpotMAX deployed 0.325 and -4.24, SpotMAX re-tuned (the
grid selection, the reference of section 7) 0.063 and -12.36. Over all nine folds: net3d 0.876, the tuned log 0.891,
SpotMAX deployed 0.694 and the grid selection 0.285; net3d and the tuned log are above the grid selection with intervals
that exclude zero and do not differ from each other. On the six registered folds the unchanged detectors reproduce
their earlier numbers exactly, while the revised net3d (0.618, from 0.878; one fold, n2rad51_HS_herm_000, overcounted
about threefold) and the revised tuned log (0.777, from 0.893; the replay now selects the configuration of the log
defaults) no longer differ from the grid selection in CCC (docs/BENCHMARK.md, "Results of the revision"). The frozen
methods change as the revised table of section 2 gives: `net3d_trainval_entry9.json` and `log_tuned_entry9.json`. The
option of this decision goes back to the project lead with these results: (1) score with the revised frozen methods, or
(2) wait for the re-export of n2rad51_HS_male_002, a further declared development revision with ten folds. Questions (b)
and (c) above stay open. No test gonad is scored until he confirms one option.

Status (2026-10-13, protocol changelog entry 16, which corrects entry 15 and the status line above): the results and
the frozen methods of the status line above are post hoc under rule 3 (protocol changelog entries 9 (j) and 15 (o)),
as is every write-up that uses them. The heat-shocked males' truth has about 17 foci per truth nucleus (17.0 to 17.9),
not per traced nucleus (there 17.1 to 19.1). The designated selection of the replayed tuned log, and 7 of its 9 fold
selections, are the configuration of the log defaults; the folds that hold out n2rad51_noHS_herm_001 and
n2rad51_noHS_male_000 select logB_d0.5_s2_k0.5_g3_m0 (z merge off), and the six-fold 0.777 is that of the fold
procedure. Three of the five designated finalists of the replay, all three at pre-gate 2, were scored on val runs that
reached the candidate cap (entry 16 (b)); the selection and the other finalist did not reach it, and the selection
stands. Under option (1) the re-export of n2rad51_HS_male_002 can no longer enter as a development revision before
test scoring: after the first scoring it could only be the revision of protocol section 10, which uses the second and
last test scoring of each method family and is post hoc. The option is still his, with questions (b) and (c) above:
no test gonad is scored until he confirms one.

Status (2026-10-13, protocol changelog entry 17): the project lead re-exported n2rad51_HS_male_002 from Imaris with
region growing off ("Should be saved"), as he said he would ("okay i will re export", entry 11). The split record and
its pinned digests now name the re-export (13,972 spots, Quality 6.847, the registered classification); its spot truth
stays out of every spot metric until a truth export of it passes the frame check, section 4.4 and the other checks of
entry 9 (d)(5), reported in a later entry, and only then do the steps that use spot truth run again with ten folds, into
folders of their own, as the further declared development revision before test scoring of entry 13 (b), post hoc under
rule 3. No test scoring has happened, so the re-export can still enter that way (entry 16 (e)). It is not an answer to
this decision: the option is still his, with questions (b) and (c) above, and no test gonad is scored until he confirms
one.

Status (2026-10-13, protocol changelog entry 18): the truth export of the re-export of n2rad51_HS_male_002 passed every
check that entry 17 declared (the frame check, the nucleus truth equal to the registered export, section 4.4 at 0.9962
against the minimum 0.9935, the checks of entry 9 (d)(5)), post hoc under rule 3. The ten-fold revision of entry 13 (b)
is prepared and declared, with the frozen methods it will give (net3d `net3d_trainval_10folds.json`; the tuned log
`log_tuned_10folds.json` if its designated selection changes, otherwise `log_tuned_entry9.json`; SpotMAX unchanged); it
runs only from the commit of entry 18, and its results go into a later entry. The project lead answered questions (c)
and (b) above as they were put to him at 12:21:01 MDT (entry 18 (b)(12); apostrophes in ASCII): "1 that's fine we have
so many datasets" to "The 2026-01-05 heat-shock males are now training data, so that session can't later serve as a test
set. OK?", so 2026-01-05 cannot be the session of a future test set (protocol section 10); and "2 sure, did you want me
to do that?" to "If someone ever re-saves the Spots of 20260105_n2_hs_1_MALE_1 in Imaris, may it join training?", which
admits that path only as put (a new save of its Spots in Imaris, then a changelog entry of its own declaring its checks,
before test scoring); the reply to him was "On 2: no, skip it for now", so none is planned and that gonad stays nucleus
truth only. On this decision he replied neither "wait" nor "score now": he re-exported n2rad51_HS_male_002 the same
morning, the course of option (2) of the entry-15 status line, and the messages of 11:57:21 and 12:21:01 MDT told him
that the one answer asked of him after the ten-fold results is "score now". Entry 18 (c)(4) (D9) declares, before any
ten-fold result: "score now" after the run scores the ten-fold frozen methods above (SpotMAX as before); the nine-fold
set of the entry-15 status line is not offered again, since the re-export has entered before any test scoring (entry 16
(e)); a choice of the nine-fold set made after the ten-fold results would be recorded verbatim and labelled as made
after them. That is a reading of his course, not his words, and it is put to him with the results. Whether and when to
score is still his: no test gonad is scored until he says "score now".

### Decision 2. Which SpotMAX configuration is the reference of section 7

Section 7: the reference is "the better of SpotMAX deployed and SpotMAX re-tuned on val (chosen on val, not
on test; both are reported on test)", under the primary metric (section 5.1 partition, Class A+B,
end-to-end). Changelog entry 2 left open which search supplies the re-tuned competitor. Package E1 scored
all 24 D3 SpotMAX sets on the two val gonads with the primary partition (trial log lines 174 to 197,
`logo_e1/val_scoring/spotmax/summary.csv`). The candidates, pooled over n2rad51_HS_herm_001 and
n2rad51_noHS_herm_001 (241 truth nuclei; with two gonads an interval only shows how far apart they are):

| candidate | parameters | val truth-mask CCC [interval] | bias | MAE | val spot F1 | how it would be chosen |
|---|---|---|---|---|---|---|
| SpotMAX deployed | triangle, radius 0.3 um, smoothing 0.08 um, effect size 3.0 | 0.556 [0.391, 0.679] | +2.91 | 2.93 | 0.428 | the 8666f07 recipe |
| D3 selection | triangle, 0.3, 0.12, 3.5 | 0.695 [0.512, 0.835] | +1.32 | 1.63 | 0.444 | the D3 rule (germline objective, three finalists on val) |
| entry 2's alternative | triangle, 0.35, 0.08, 4.0 | 0.723 [0.652, 0.762] | +1.00 | 1.68 | 0.420 | the best D3 finalist on val by the harness's older traced subset |
| best D3 set on val under the primary metric | otsu, 0.3, 0.08, 4.0 | 0.822 [0.734, 0.822] | -0.12 | 1.17 | 0.365 | section 6.4 read literally over the D3 search space (every set has a val score now) |

Every re-tuned candidate beats the deployed recipe on val under the primary metric, so by section 7 the
reference is a re-tuned configuration; which one is the decision. Points to weigh: the otsu set is the
literal 6.4 choice within the D3 space but has the lowest F1 of the four, in D3 every otsu set starved the
heat-shock training gonad (bias -0.88 to -1.41 foci per germline nucleus on n2rad51_HS_herm_000), and held out
in the leave-one-gonad-out estimate that gonad got a truth-mask bias of -2.33 from it; the leave-one-gonad-out
estimate of "SpotMAX re-tuned" measures the D3 procedure, which picked the D3 selection in four of six folds and
the otsu set in the two heat-shock folds; none of the re-tuned candidates was chosen under the primary metric by
a rule written before the search. Parameter files: `final_params/spotmax_retuned_primary.json` (sha256 84e27f45...) for the otsu set,
`spotmax_retuned_d3.json` for the D3 selection. In the leave-one-gonad-out estimate (docs/BENCHMARK.md,
"Leave-one-gonad-out") the re-tuning procedure (the D3 search replayed per fold) reached a truth-mask CCC of
0.742 [0.655, 0.812] against 0.601 [0.499, 0.718] for the deployed recipe; that difference, +0.141 [-0.043,
+0.292], does not exclude zero.

The section 8 item 2 grid (thresholding method in otsu, li, triangle, yen; `effect_size_min` 2.0 to 4.0;
`z_merge_gap_um` 0.6, 0.8, 1.0; 60 configurations) has never been run. Running it now is allowed (development
on train and val, val score for every configuration, section 6.4), takes several hours of SpotMAX, and would
replace the D3 search as the re-tuned competitor if the project lead prefers the pre-registered grid.

Status (2026-10-11, protocol changelog entry 6): a choice, which the reply "Yes to everything and merged" does not
settle by itself. Entry 6 declares the section 8 item 2 grid before it is run: 60 sets
(`pg_<method>_e<effect size>_z<gap>`, search `prereg_spotmax_grid`) scored on the four usable train gonads and the
two val gonads, selected on val by the pooled truth-mask CCC rounded to 3 decimals, then the spot F1 inside R,
then listing order. Recommendation: that selection becomes SpotMAX re-tuned, and by section 7 the reference, once
the project lead confirms it and the readings of entry 6 (c) and (d); the three D3 rows above stay reported as post
hoc. The grid contains the otsu row above (`pg_otsu_e4.0_z0.8`, in effect the same set), so if its val score
reproduces, the selection is that set unless a grid set not scored yet ranks above it on val. The choice is made on
train and val only (section 5); entry 6 (b) lists what was known before, including a pre-protocol sweep on two
test gonads.

Status (2026-10-11, protocol changelog entry 7): still the project lead's choice. The grid ran as entry 6 declared
it, on train and val only, and every check of entry 6 (c) passed. Its selection is `pg_otsu_e4.0_z0.8`
(`threshold_otsu`, effect size 4.0, z gap 0.8), in effect the otsu row of the table above: val truth-mask CCC 0.822
against 0.556 for the deployed recipe (rank 23 of 60). By section 7 it is the reference once he confirms it and the
readings of entry 6 (c) and (d); the readings that can be recomputed (val only, rounding, early stop) leave the set
unchanged. Its scores per gonad show the heat-shock concern above: truth-mask bias -2.33 on the train gonad
n2rad51_HS_herm_000 (CCC 0.587; the deployed recipe +0.91, CCC 0.843) and -0.77 on the val gonad
n2rad51_HS_herm_001, against +0.72 to +1.38 on the three train gonads without heat shock and +0.30 on the val gonad
n2rad51_noHS_herm_001. Computed after the run: every grid set misses one of the two heat-shock hermaphrodites by at
least 1.7 foci per truth nucleus. Two of the four test gonads are heat-shocked, and no heat-shock male informed
this choice. If the selection undercounts them as it undercounts both heat-shock hermaphrodites here, it is a weak
reference on half the test set, which makes W1 and W2 easier to meet. In the leave-one-gonad-out estimate the grid
procedure reaches a truth-mask CCC of 0.684 [0.551, 0.788], against 0.601 for the deployed recipe and 0.742 for the
D3 procedure, and neither difference excludes zero. The D3 procedure is ahead on the F1 inside R (+0.024 [+0.002,
+0.041]) and on the detector-only CCC (+0.041 [+0.011, +0.074]) and mean absolute error (-0.210 [-0.348, -0.088]),
intervals that exclude zero (docs/BENCHMARK.md, "Results of the pre-registered SpotMAX grid"). The recommendation
of entry 6 is kept: section 8 item 2 defines SpotMAX re-tuned as this grid chosen by the val rule, the grid ran and
selected as declared, and the selection is the set entry 6 expected, with the heat-shock undercount this decision
already named. Another set chosen now would be chosen after these results. If he picks another set, the pre-scoring
entry records it and its reason as post hoc.

Status (2026-10-11, protocol changelog entry 8): settled. The project lead's reply of 2026-10-11 at 10:01 MDT
(16:01 UTC), which begins "Grid sure", answers the request to reply "grid" for this decision or to name another
setting. The grid selection `pg_otsu_e4.0_z0.8` is frozen as SpotMAX re-tuned
(`final_params/spotmax_retuned_grid.json`, section 2), and by section 7 it is the reference: on val its truth-mask
CCC is 0.822 against 0.556 for the deployed recipe, which is one of the 60 sets, so the selection cannot rank below
it there. That comparison is in sample, since the set was chosen on the same two gonads. The readings of
entry 6 (c) and (d), which entry 7 (g) and the description of PR #10 put to him in writing with the recommendation,
are read as confirmed with it; that is a reading of the two words "Grid sure", not his words. The three D3 rows of
the table above stay reported as post hoc and are not candidates. The otsu row among them runs the same
configuration as the grid selection; what is frozen is the grid's own file. The heat-shock concern above stands:
the reference undercounts both heat-shock hermaphrodites of train and val, so the write-up reports that caveat next
to W1 and W2, which are evaluated exactly as section 7 states them (entry 8 (d) gives the caveat in full). Not
settled by the reply: protocol section 10 says that if W1 fails "the SpotMAX reference stays deployed", but the
pipeline runs the 8666f07 recipe, not this set. Which setting it keeps in that case is for the project lead to
settle before test is scored, and the pre-scoring entry records it (section 4.2).

### Decision 3. The limit of two test scorings per method family

`scripts/bench_detector.py --final` refuses a third test scoring of the same backend (protocol section 10:
"at most twice per method family"). SpotMAX deployed and SpotMAX re-tuned are both backend `spotmax`, so
scoring both (section 7 asks for both) uses both scorings of that family and leaves no declared revision for
either. Either accept that (the SpotMAX family is the reference, not a method under development), or decide
before scoring that they are two families and give the guard a family label (a code change, declared in the
pre-scoring changelog entry). This document assumes the first.

Status (2026-10-11, protocol changelog entry 6): accepted as written by the reply "Yes to everything and merged"
(the question put it as "accept that"); SpotMAX deployed and SpotMAX re-tuned use the two test scorings of the
`spotmax` family.

### Decision 4. The competitors of section 8 that do not exist yet

Spotiflow (zero-shot and fine-tuned), RS-FISH, Imaris with one global Quality threshold and the human blind
subset (section 9) are not implemented. The win conditions compare a new detector with the SpotMAX reference
only, so they do not block W1 and W2; but each later competitor is then scored after the first results are
known, which the write-up must say.

Status (2026-10-11, protocol changelog entry 6): accepted as written; these competitors are scored after the first
test results, and the write-up says so.

### Decision 5. The post hoc choices that stand

Changelog entry 2 (the D3 searches: germline objective, finalists only on val, the SpotMAX grid not run) and
entry 3 (net3d's selection rule) are post hoc under protocol rule 3, and so is the leave-one-gonad-out design
of entry 4. The write-up labels them. Confirm that the methods below are scored as they stand.

Status (2026-10-11, protocol changelog entry 6): accepted with a reading. The choices of entries 2, 3 and 4
stand and are labelled post hoc, and the methods below are scored as they stand except SpotMAX re-tuned:
decision 2 leaves that row open, and "the SpotMAX grid not run" no longer holds once entry 6 runs the grid.
The SpotMAX re-tuned row of section 2 stands only if decision 2 keeps the D3 selection.

Status (2026-10-11, protocol changelog entry 8): decision 2 is settled and did not keep the D3 selection. The
SpotMAX re-tuned row of section 2 is now the grid selection, and the grid ran (entry 7), so decision 2 leaves no
row of section 2 open. Decision 1 is still open: if the re-exported heat-shock males enter training, the methods it
names change in a declared revision before any test scoring. The two D3 parameter files stay in `final_params/` for
the record and are not scored.

Status (2026-10-12, protocol changelog entry 9): the net3d row, and the tuned log row if its replayed selection
differs, are replaced by entry-9 files before test scoring; both are post hoc. The SpotMAX rows stand.

Status (2026-10-13, protocol changelog entry 15): both rows are replaced. The replayed D3 procedure selects
logA_d0.5_s2 (in effect the log defaults), not the set of `log_tuned.json`, so `log_tuned_entry9.json` replaces it, and
`net3d_trainval_entry9.json` (the model trained on nine gonads) replaces `net3d_trainval.json`; both are post hoc
(entries 9 and 15). The SpotMAX rows stand.

Status (2026-10-13, protocol changelog entry 16): the tuned log of the revised table is post hoc under entries 9, 11, 12
and 15, like net3d, since its replay tuned on n2rad51_HS_male_000 and n2rad51_20260105a_HS_male_001; the status line
above cited entries 9 and 15 only.

## 2. The frozen methods

| competitor | backend | parameters (`--params` file) | weights | what it is |
|---|---|---|---|---|
| SpotMAX deployed | spotmax | `final_params/spotmax_deployed.json`: `{}` (sha256 e3566b3a...) | none | the 8666f07 recipe: the golden runs' `spots:` section, config hash 84bdda5a5230 |
| SpotMAX re-tuned (grid selection `pg_otsu_e4.0_z0.8`; decision 2) | spotmax | `final_params/spotmax_retuned_grid.json` (sha256 919c4ab5...) | none | `threshold_otsu`, effect size 4.0, z merge gap 0.8 um, every other `spots:` key as deployed (radius 0.3 um, smoothing 0.08 um): the section 8 item 2 grid, selected on val by the rule of protocol changelog entry 6 (d), whose readings are post hoc (entries 6 to 8), and, by section 7, the SpotMAX reference |
| tuned log | log | `final_params/log_tuned.json` (sha256 3a65e844...) | none | diameter 0.5 um, SNR 2.5, pre-gate 2, z merge off (profile `bench_log_tuned`) |
| net3d (train plus val) | net3d | `final_params/net3d_trainval.json` (sha256 2f9097e6...) | `models/net3d_trainval_20261010/model.pt`, sha256 d091b9ef58850783387a4f5aa0fea5ecc9310abc027145771936fecd0de001c0 (commit 240cc72) | the frozen recipe of net3d_20261010 (epoch 4 of the 40-epoch schedule, threshold 0.40) retrained on the six train and val gonads with usable spot truth (protocol section 6.1) |

Entry 9: the net3d and tuned log rows hold until the entry-9 revision runs; `net3d_trainval_entry9.json` (a new
model sha256) and, if written, `log_tuned_entry9.json` replace them, and the pre-scoring entry names them. The
revision does not touch the SpotMAX rows.

`final_params/` is `C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/final_params/`; its six files have
the sha256 given here. Five were written on 2026-10-10: the three files above other than SpotMAX re-tuned's, and
two D3 files that stay for the record and are not scored, `spotmax_retuned_d3.json` (the D3 selection, sha256
67db0cf7...) and `spotmax_retuned_primary.json` (the D3 otsu set, sha256 84e27f45...). `spotmax_retuned_grid.json`
was written on 2026-10-11, after decision 2 was settled; protocol changelog entry 8 gives its full sha256. SpotMAX
re-tuned is frozen as the grid selection `pg_otsu_e4.0_z0.8` (entries 7 and 8), and the file holds exactly its
three keys, the `selected_params` under `designated` in `prereg_grid/selection.json`:
`{"thresholding_method": "threshold_otsu", "effect_size_min": 4.0, "z_merge_gap_um": 0.8}`. Laid over the golden
`spots:` section (`detector.effective_section`), these keys give the same section as
`spotmax_retuned_primary.json`: the two files run the same SpotMAX configuration under different parameter hashes
(dc085381... for the grid file, the hash of this set in the grid's trial records; 35d0c2a3... for
`spotmax_retuned_primary.json`). The grid's own file is the one scored, because entry 6 (d) makes the D3 sets
non-candidates now that the grid is confirmed.

Revised after the entry-9 revision (protocol changelog entry 15, 2026-10-13). The table above is kept as it was
written; for the test scoring its tuned log and net3d rows are replaced by those below.

| competitor | backend | parameters (`--params` file) | weights | what it is |
|---|---|---|---|---|
| SpotMAX deployed | spotmax | `final_params/spotmax_deployed.json`: `{}` (sha256 e3566b3a..., unchanged) | none | the 8666f07 recipe, config hash 84bdda5a5230; not touched by the revision |
| SpotMAX re-tuned (grid selection `pg_otsu_e4.0_z0.8`; decision 2) | spotmax | `final_params/spotmax_retuned_grid.json` (sha256 919c4ab5..., unchanged) | none | selected on the two val gonads (protocol changelog entry 6 (d)); the additions are train gonads and enter no val selection (entry 9 (c)), so it does not change; by section 7 the SpotMAX reference |
| tuned log | log | `final_params/log_tuned_entry9.json` (sha256 0fa43134d1e66e5221755e99d5af461d7e80ab55e6cbd747242bebad5c5338eb): `{"log": {"diameter_um": 0.5, "min_snr": 2.0}}` | none | logA_d0.5_s2, the designated selection of the D3 procedure replayed with n2rad51_HS_male_000, n2rad51_20260105a_HS_male_001 and _002 in train (entry 9 (g)): diameter 0.5 um and SNR 2.0, every other key at its default (pre-gate 3.0; z merge as in the golden `spots:` section), so in effect the log defaults; post hoc (entries 9, 11, 12 and 15). Replaces `log_tuned.json` (3a65e844...), which stays in `final_params/` and is not scored |
| net3d (train plus val) | net3d | `final_params/net3d_trainval_entry9.json` (sha256 953b75193cbaee28d60460bdf9f4d81e87c4c8808991563116716eca064ea533) | `models/net3d_trainval_entry9/model.pt`, sha256 ef986d09671907918cb485e7a13363daf25eba0c051b7456d6fcd5032a0eef52 (commit 08f6285) | the frozen recipe of net3d_20261010 (epoch 4 of the 40-epoch schedule, threshold 0.40, seed 0) retrained on the nine train and val gonads with usable spot truth, in this order: n2rad51_HS_herm_000, n2rad51_noHS_herm_000, n2rad51_noHS_male_000, n2rad51_noHS_male_001, n2rad51_HS_herm_001, n2rad51_noHS_herm_001, n2rad51_HS_male_000 (by the recorded export of entry 11), n2rad51_20260105a_HS_male_001 and _002; post hoc (entries 9, 11, 12 and 15). Replaces `net3d_trainval.json` (2f9097e6..., model d091b9ef...), which stays and is not scored |

`final_params/` now holds eight files: the six of 2026-10-10 and 2026-10-11, unchanged, and the two entry-9 files,
written by the launcher of the revision on 2026-10-12 (`log_tuned_entry9.json` at 22:53 MDT, `net3d_trainval_entry9.json`
at 23:51 MDT).

## 3. The order

1. Pre-flight (no metric is computed).
2. The pre-scoring changelog entry (protocol section 10: code commit, config hash and weight sha256 are
   written before test is scored), committed.
3. SpotMAX deployed, then the golden check of its detections.
4. SpotMAX re-tuned: the grid selection `pg_otsu_e4.0_z0.8` (decision 2, settled in protocol changelog entry 8;
   `final_params/spotmax_retuned_grid.json`), by section 7 the SpotMAX reference.
5. Tuned log.
6. net3d (train plus val), on the GPU under the workstation's lock.
7. The comparison (pooling only) and the win conditions.
8. The post-scoring changelog entry, committed.

Every method is frozen before the first is scored, so no test number can shape another method; the two SpotMAX
configurations are scored first only because the reference, SpotMAX re-tuned (decision 2), is one of them and
everything else is compared with it.

## 4. The commands

Run from the repository root of the commit named in the pre-scoring entry, in PowerShell with the
repository's environment active (`python` is its interpreter and imports this checkout's `germquant`).
That holds as written only in the main checkout: the environment's editable install points at the main checkout's
src, so from any other checkout (a worktree) set `$env:PYTHONPATH = "<that checkout>/src"` first; then the
`germquant` executable, `python -m germquant` and `python -c "import germquant ..."` all import that checkout
(protocol changelog entries 13 (c)(1) and 14 (g)). The pre-flight prints where the package comes from.

### 4.1 Pre-flight

```powershell
git status --porcelain                 # must print nothing (--final refuses uncommitted code)
git rev-parse HEAD                     # goes into the pre-scoring entry
# germquant must come from this checkout's src (not the main checkout's, unless this is the main checkout)
python -c "import germquant; print(germquant.__file__)"
$env:CUDA_VISIBLE_DEVICES = "-1"
# the one test gonad without a bench cache: images and labels only, no metric
python scripts/bench_cache.py --tags n2rad51_noHS_male_002 --verify-sha256
# the golden run of that gonad must be the deployed pipeline (git sha 8666f07..., config hash 84bdda5a5230)
python -c "import json; m = json.load(open(r'C:/Users/ryane/germquant_golden/8666f07/n2rad51_noHS_male_002/run_manifest.json')); print(m['git_sha'], m['config_hash'])"
# the net3d weights are the frozen ones: since the entry-9 revision the model of net3d_trainval_entry9.json, which must
# print ef986d09671907918cb485e7a13363daf25eba0c051b7456d6fcd5032a0eef52 (protocol changelog entry 15; the earlier
# net3d_trainval_20261010 model, d091b9ef..., is not scored)
python -c "from germquant.provenance import file_sha256; print(file_sha256(r'C:/Users/ryane/germquant_golden/models/net3d_trainval_entry9/model.pt'))"
# the additions of protocol changelog entry 9 are the recorded files (must print [])
python -c "from germquant.validate import benchmark_split as bs; print(bs.check_addition_files(bs.load_additions('tests/golden/rad51_additions.json'), hashes=True))"
# the recorded spot source of protocol changelog entry 11 (n2rad51_HS_male_000's 2025-12-31 export, on which the
# train-plus-val net3d and the tuned log's replayed selection depend) is the pinned file, checked against the split
# (must print [])
python -c "from germquant.validate import benchmark_split as bs; print(bs.check_spot_source_files(bs.load_spot_sources('tests/golden/rad51_spot_sources.json', bs.load_split('tests/golden/rad51_split.json')), hashes=True))"
# the parameter files are the frozen ones (sha256 as in section 2; Get-FileHash prints them in upper case):
# spotmax_retuned_grid.json, the grid selection of decision 2, is 919C4AB5... (protocol changelog entry 8);
# log_tuned_entry9.json is 0FA43134... and net3d_trainval_entry9.json 953B7519... (protocol changelog entry 15)
Get-FileHash -Algorithm SHA256 C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/final_params/*.json
# protocol section 4.4 on the test gonads: the geometric truth assignment against the lab's mask channel;
# agreement should match the train and val gonads (99.35 to 99.99 percent, read as a minimum of 99.35 by
# protocol changelog entry 10 (b)); there the differences were boundary spots, spots whose mask-channel value is a
# removed mask fragment (the larger group), and 3 spots that one region-growing gonad assigns to different nuclei
python scripts/check_truth_assignment.py --tags n2rad51_HS_herm_002 n2rad51_HS_male_001 n2rad51_noHS_herm_002 n2rad51_noHS_male_002
# the truth frame of each test gonad (protocol changelog entry 10 (e)): Imaris's per-spot Intensity Center test
# of entry 9 (d)(3) only, no metric; the export goes to a scratch folder and is not used. Every line must print
# "stored" (the frame the test gonads' truth was exported with); any other answer stops the scoring until that
# gonad's truth is exported again with the chosen frame, through a changelog entry
$S = Get-Content tests/golden/rad51_split.json -Raw | ConvertFrom-Json
foreach ($g in ($S.gonads | Where-Object { $_.role -eq "test" })) { germquant imaris-truth $g.ims.path --nd2 $g.nd2.path --frame-check --pre-flight --no-masks --out "$env:TEMP/preflight_frame/$($g.tag)" }
python -c "import glob, json; [print(p, (json.load(open(p))['frame_check'].get('chosen') or {}).get('name')) for p in sorted(glob.glob(r'$env:TEMP/preflight_frame/*/*__imaris_params.json'))]"
# no test scoring has happened
Test-Path C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/final_runs.jsonl   # must be False
```

The golden run of n2rad51_noHS_male_002 was made by the pinned commit 8666f07 (finished 2026-10-10 after 205.5
minutes in the flooded spot stage), so the deployed competitor on that gonad is the 8666f07 pipeline; the
check above confirms it. The same gonad rerun by the code with the KD-tree peak filter is byte-identical to it.

### 4.2 The pre-scoring changelog entry

Append to section 11 of docs/BENCHMARK_PROTOCOL.md (fill the brackets), then commit:

> N. [date]. Test scoring declared before any test result (protocol section 10). Methods, frozen at commit
> [git sha]: SpotMAX deployed (backend spotmax, params `{}`, the golden `spots:` section, config hash
> 84bdda5a5230); SpotMAX re-tuned (backend spotmax, the grid selection `pg_otsu_e4.0_z0.8`, params file
> `spotmax_retuned_grid.json`, sha256 919c4ab5...); the SpotMAX reference of section 7 is SpotMAX re-tuned, chosen
> on val under the primary metric (truth-mask CCC 0.822 against 0.556 for the deployed recipe, a comparison in
> sample; decision 2, entry 8); if W1 fails, SpotMAX stays deployed with [the 8666f07 recipe or the grid selection]
> (section 10); tuned log (backend log, params file `[log_tuned_entry9.json]`, sha256 [0fa43134...]:
> [logA_d0.5_s2, in effect the log defaults], the designated selection of the D3 procedure replayed with the
> heat-shocked males in train; config hash 84bdda5a5230 (the golden runs' config, whose `spots:` section its params
> are laid over), parameter hash [as its run records give it]; post hoc, protocol changelog entries 9, 11, 12 and
> 15[, and the entry of a later revision]); net3d (backend net3d, params file `[net3d_trainval_entry9.json]`, sha256
> [953b7519...], naming `[models/net3d_trainval_entry9/model.pt]`, sha256 [ef986d09...], trained at [08f6285];
> config hash 84bdda5a5230 (as for the tuned log), parameter hash [as its run records give it]; the frozen recipe of
> net3d_20261010 retrained on the [nine or ten] train and val gonads with usable spot truth ([seven or eight]
> registered: the six with spot truth from their .ims, n2rad51_HS_male_000 by the recorded export of entry 11
> (`tests/golden/rad51_spot_sources.json`, sha256 [8e145fb1...]; its criteria post hoc, entry 12 (a))[,
> n2rad51_HS_male_002 by its re-export of entry N]; [two] additions of entry 9, as the decisions file of entries 10
> and 11 and any later entry keeps them); post hoc, protocol changelog entries 9, 11, 12 and 15[, and the entry of a
> later revision]).
> Decisions 1 to 5 of docs/FINAL_SCORING.md: [one line
> each]. No test result has been computed at this point.

### 4.3 SpotMAX deployed

```powershell
$env:CUDA_VISIBLE_DEVICES = "-1"
$F = "C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07"
python scripts/bench_detector.py --backend spotmax --params $F/final_params/spotmax_deployed.json `
    --split test --final --detector-only --jobs 2 `
    --reason "one-time test scoring (protocol section 10): SpotMAX deployed, the 8666f07 recipe" `
    --out $F/final/spotmax_deployed
```

Two SpotMAX workers need about 25 GB of memory; n2rad51_noHS_male_002 may take much longer than the others
(its detector flooded before the peak-filter fix). Then check, with no metric, that the detections are the
golden spot tables (protocol section 8 item 1: the deployed competitor's outputs are the goldens):

```powershell
@'
import pandas as pd
from germquant.benchmark import detector as bd
from germquant.benchmark.manifest import load_manifest
E = {e.tag: e for e in load_manifest(bd.DEFAULT_MANIFEST)}
for t in ["n2rad51_HS_herm_002", "n2rad51_HS_male_001", "n2rad51_noHS_herm_002", "n2rad51_noHS_male_002"]:
    run = pd.read_csv(f"C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/final/spotmax_deployed/detections/{t}__spots.csv", float_precision="round_trip")
    gold = bd.read_golden_spots(E[t].results_dir, bd._entry_stem(E[t]))
    print(t, bd.compare_with_golden(run, pd.DataFrame(), gold)["identical"])
'@ | python -
```

Every line must say True. A False means the scored detections are not the goldens; stop and record it.

### 4.4 SpotMAX re-tuned

```powershell
python scripts/bench_detector.py --backend spotmax --params $F/final_params/spotmax_retuned_grid.json `
    --split test --final --detector-only --jobs 2 `
    --reason "one-time test scoring (protocol section 10): SpotMAX re-tuned, grid selection pg_otsu_e4.0_z0.8" `
    --out $F/final/spotmax_retuned
```

This is the second and last `--final` run the guard allows for backend `spotmax` (decision 3). The `--params` file
is the grid selection's own (decision 2, settled in protocol changelog entry 8; sha256 919c4ab5..., section 2), and
the pre-scoring entry names it. Its run records will carry the parameter hash dc085381..., the hash of this set in
the grid's trial records. `spotmax_retuned_primary.json` runs the same configuration under another hash and is not
used.

### 4.5 Tuned log

```powershell
python scripts/bench_detector.py --backend log --params $F/final_params/log_tuned_entry9.json `
    --split test --final --detector-only --jobs 2 `
    --reason "one-time test scoring (protocol section 10): tuned log (the D3 procedure replayed after the entry-9 revision)" `
    --out $F/final/log_tuned
```

Since the entry-9 revision the `--params` file is `log_tuned_entry9.json` (sha256 0fa43134..., section 2, revised
table; protocol changelog entry 15); `log_tuned.json` stays in `final_params/` and is not scored.

### 4.6 net3d (train plus val)

```powershell
Remove-Item Env:CUDA_VISIBLE_DEVICES -ErrorAction SilentlyContinue
$env:PYTORCH_CUDA_ALLOC_CONF = "expandable_segments:True"     # a no-op on Windows, kept for the workstation rule
python C:/Users/ryane/germquant_golden/gpu_lock.py run --name "final test scoring: net3d_trainval" -- `
    python scripts/bench_detector.py --backend net3d --params $F/final_params/net3d_trainval_entry9.json `
    --split test --final --detector-only `
    --reason "one-time test scoring (protocol section 10): net3d retrained on the nine train and val gonads by the frozen recipe" `
    --out $F/final/net3d_trainval
```

Since the entry-9 revision the `--params` file is `net3d_trainval_entry9.json` (sha256 953b7519...), which names
`models/net3d_trainval_entry9/model.pt` (sha256 ef986d09..., section 2, revised table; protocol changelog entry 15);
`net3d_trainval.json` and its model stay and are not scored. The output folder and the report label keep the name
`net3d_trainval`.

### 4.7 The comparison

Pooling only; every number comes from the four runs above. New methods first, so that every paired
difference reads new minus reference:

```powershell
$env:CUDA_VISIBLE_DEVICES = "-1"
python scripts/bench_logo.py report `
    --tags n2rad51_HS_herm_002 n2rad51_HS_male_001 n2rad51_noHS_herm_002 n2rad51_noHS_male_002 `
    --heldout net3d_trainval=$F/final/net3d_trainval log_tuned=$F/final/log_tuned `
              spotmax_retuned=$F/final/spotmax_retuned spotmax_deployed=$F/final/spotmax_deployed `
    --out $F/final/report
```

`report` refuses a folder that lacks a gonad. A method that produced no output for a gonad (its run lists the
gonad under `skipped`) gets a row saying so in every table (protocol section 10); the win conditions are then
evaluated over the gonads where both methods produced output: rerun `report` with `--tags` holding those
gonads for that pair and label the result with their number (3 of 4, say).

Reading the win conditions (protocol section 7) from `final/report/logo_paired.csv`, with R the reference of
decision 2, SpotMAX re-tuned (`spotmax_retuned` in the report):

* W1 (counts): the row with `pair` = "net3d_trainval - R", `block` = truth_masks, `stat` = ccc: `diff` at least
  0.03 and `lo` above 0 (the paired cluster bootstrap by gonad, 1000 resamples, seed 20251105). The same row
  for "log_tuned - R" says whether the tuned log wins.
* W2 (spots): `block` = spots_truth_masks (protocol section 5.2: matching restricted to R, the truth masks;
  spots_all is the whole-field score, kept for comparison only), `stat` = f1: `diff` not below -0.02 (point estimates; the interval is
  reported). `logo.DETECTION_BLOCKS` holds this block since commit ca8537d (`"spots_truth_masks": "tm_spots_"`,
  with `spots_truth_masks_class_a` and `spots_truth_masks_detector_only`), so `bench_logo.py report` pools and
  pairs it; a held-out folder without the `tm_spots_` columns (written before commit e19a029) is left out of these
  blocks (protocol changelog entry 7). If `logo_paired.csv` has no spots_truth_masks row for a pair, W2 cannot be
  read for it: stop and record it.
* W3 does not apply: no method changes the nucleus segmentation.
* Fragility (section 6.2): for every pair, compare the direction of the test difference with the
  leave-one-gonad-out difference of docs/BENCHMARK.md ("Leave-one-gonad-out"); a test result whose direction
  the leave-one-gonad-out estimate contradicts is reported as fragile. Once the entry-9 revision has run, the
  directions are read from its leave-one-gonad-out report (all folds, with the six registered folds and the
  heat-shocked males reported next to it; docs/BENCHMARK.md, "Entry 9"), not from the six-fold report below. The leave-one-gonad-out directions, for
  reference: net3d and the tuned log are above both SpotMAX configurations on the truth-mask CCC, the germline
  and traced CCC and spot F1 (the paired intervals exclude zero); net3d against the tuned log shows no
  difference (truth-mask CCC -0.015 [-0.084, +0.073]); SpotMAX re-tuned against deployed shows none except a
  smaller mean absolute error. These SpotMAX re-tuned directions are those of the D3 replay. Decision 2 made the
  grid selection the re-tuned competitor (protocol changelog entry 8), so its pairs are read from the grid replay
  of entries 6 and 7 instead. Both replays, with the section 5.2 F1 inside R for every pair, are in
  `prereg_grid/logo/report_f1R/logo_paired.csv` (docs/BENCHMARK.md, "Results of the pre-registered SpotMAX grid").
  That file's labels are not the test report's: there the grid replay is `spotmax_grid`, the D3 replay is
  `spotmax_retuned` and net3d is `net3d`. So the test pair "net3d_trainval - spotmax_retuned" is read against
  "net3d - spotmax_grid", "log_tuned - spotmax_retuned" against "log_tuned - spotmax_grid", and
  "spotmax_retuned - spotmax_deployed" against "spotmax_deployed - spotmax_grid" with its sign reversed. In that
  report net3d and the tuned log are above all three SpotMAX configurations on the F1 inside R (+0.221 to +0.255,
  the intervals exclude zero) and better than the grid procedure on every CCC, the truth-mask mean absolute error
  (lower) and both spot F1. On the F1 inside R no difference survives between net3d and the tuned log (+0.010
  [-0.058, +0.094]) or between the D3 replay and the deployed recipe (+0.014 [-0.073, +0.069]). The grid procedure
  minus the deployed recipe: truth-mask CCC +0.083 [-0.079, +0.197], F1 inside R -0.011 [-0.085, +0.046], mean
  absolute error -1.016 [-1.518, -0.344], whole-field F1 -0.049 [-0.098, -0.005].
* Report next to the primary numbers: Class A and Class B truth (truth_masks_class_a; per_gonad.csv
  `tm_class_b_*`), the detector-only count (truth_masks_detector_only), the germline and traced partitions,
  and spot F1 without Class B, every one per gonad and pooled (section 5.4), with the run records of
  final_runs.jsonl (versions, git sha, params sha256) and the split file sha256.

The fragility directions after the entry-9 revision (protocol changelog entry 15), post hoc under rule 3 (entries 9 (j),
15 (o) and 16 (d)), as the fragility item above asks: they are read from `entry9/logo/report/logo_paired.csv` (all nine
folds; `report_registered6` and `report_hs_males` are reported next to it). Its labels are net3d (the fold models),
log_tuned, log_defaults, spotmax_deployed and spotmax_grid (the grid selection; it has no D3 SpotMAX replay). Each test
pair is read against the pair below, with the truth-mask CCC and the F1 inside R of that report (an asterisk marks an
interval that excludes zero):

| test pair | read against (entry-9 report) | truth masks CCC | F1 inside R |
|---|---|---|---|
| net3d_trainval - spotmax_retuned | net3d - spotmax_grid | +0.591 [+0.083, +0.807] * | +0.419 [+0.218, +0.546] * |
| log_tuned - spotmax_retuned | log_tuned - spotmax_grid | +0.606 [+0.156, +0.776] * | +0.375 [+0.203, +0.488] * |
| net3d_trainval - log_tuned | net3d - log_tuned | -0.015 [-0.181, +0.096] | +0.044 [-0.071, +0.111] |
| spotmax_retuned - spotmax_deployed | spotmax_deployed - spotmax_grid, sign reversed | -0.409 [-0.589, -0.018] * | -0.079 [-0.131, -0.007] * |
| net3d_trainval - spotmax_deployed | net3d - spotmax_deployed | +0.182 [-0.021, +0.438] | +0.340 [+0.184, +0.442] * |
| log_tuned - spotmax_deployed | log_tuned - spotmax_deployed | +0.198 [+0.099, +0.391] * | +0.296 [+0.180, +0.370] * |

On the six registered folds the truth-mask CCC of net3d and of the tuned log no longer differs from the grid selection
(-0.066 [-0.258, +0.277] and +0.093 [-0.117, +0.342]); the write-up reports that next to the all-fold directions.

### 4.8 The post-scoring changelog entry

Append to section 11 (fill the brackets), then commit:

> N+1. [date]. Test scored once per method (final_runs.jsonl lines [n to m], commit [git sha]): SpotMAX
> deployed, SpotMAX re-tuned, tuned log, net3d (train plus val), on [4 of 4] test gonads [list any gonad
> without output]. Primary metric (section 5.1 partition, Class A+B, end-to-end), pooled CCC [95 percent
> interval]: SpotMAX deployed [..], SpotMAX re-tuned [..], tuned log [..], net3d [..]. W1 net3d against the
> reference R (SpotMAX re-tuned, decision 2): difference [..] [interval]: [met / not met]; W2: F1 difference [..]:
> [met / not met]. Tuned log against R: W1 [..], W2 [..]. Decision (section 7): [the new detector replaces SpotMAX
> / SpotMAX stays deployed (section 10) with [the setting the pre-scoring entry names] / mixed result, not
> replaced]. Fragile results (section 6.2): [none / which].
> These are the first test results of the new methods; any change from here on is post hoc (rule 3) and a
> second scoring needs a revision declared here first (section 10).

## 5. What must not happen

* No method is changed, tuned or re-selected after any test number is seen; a revision is declared in the
  changelog before its (second and last) scoring.
* No test gonad is scored with `--no-trial-log`, with `scripts/benchmark_spots.py`, or outside `--final`.
* The SpotMAX reference is not chosen, or changed, on test numbers.
