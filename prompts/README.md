# LLM prompts

Every prompt used to build and evaluate CrossValQA / MutQA, verbatim, with the
model and decoding parameters each was run with. Reproduced here as a single
reference. Prompts 1–3 are also in the pipeline source files named beneath
them; prompts 4–6 come from the paper's evaluation scripts, which are not part
of this repository, so their file names are given for provenance only.

Six roles, in pipeline order:

| # | Role | Model | Output |
|:--|:--|:--|:--|
| 1 | QA generation | `google/gemini-2.5-flash` | 2–3 QA pairs + rationale + evidence IDs |
| 2 | Independent answering | `deepseek/deepseek-chat` | answer + rationale + confidence |
| 3 | Cross-validation judge | `meta-llama/llama-3.1-8b-instruct` | binary `supported` |
| 4 | QA-pair auditing judge | `google/gemini-2.0-flash-001`, `anthropic/claude-3.5-haiku` | 0–2 correctness |
| 5 | MutaPLM evaluation judge | `deepseek/deepseek-v3.2`, `google/gemini-2.0-flash-001`, `openai/gpt-4o-mini` | 0–5 × 3 dimensions |
| 6 | Retrieval-experiment judge | `google/gemini-2.5-flash-lite`, `anthropic/claude-haiku-4.5`, `deepseek/deepseek-v3.2` | 0/1/2 + abstention |

Roles 1–3 build the dataset; 4–6 evaluate it. Information isolation is the point
of the design: role 2 never sees role 1's answer or rationale, and role 3 sees
only (answer, rationale) pairs without knowing which system produced them.

---

## 1. QA generation (Phase 1)

`claim_extractor/phase_1/phase1_generator_qa.py` · `google/gemini-2.5-flash`

**Inputs:** mutation identifier, evidence sentences extracted from the article.
**Output:** JSON list of QA pairs, each with `question`, `answer`, `rationale`,
`evidence_sentence_ids`, `question_type`.

```
SYSTEM:
You are a scientific QA generator specializing in genetics and molecular biology.

Your task is to generate high-quality question-answer pairs about genetic mutations
based on evidence from research papers.

CRITICAL RULES:
1. ALWAYS include the mutation identifier (e.g., c.125C>T, p.Thr42Met) in BOTH the
   question AND answer
2. Generate 2-3 diverse questions covering different aspects (pathogenicity, disease
   association, functional effect, mechanism, clinical features)
3. Each answer must be directly supported by the provided evidence sentences
4. Each rationale must:
   - Quote at least one relevant phrase from the evidence (in quotation marks)
   - Cite the specific sentence number(s) using [N] format
   - Explain HOW the quote supports the answer
5. If the evidence supports a negative finding (e.g., "does NOT cause...", "is NOT
   associated with..."), include it as one of the questions. Do not force negative
   framing if unsupported by evidence.
6. Do NOT infer information not explicitly stated in the evidence
7. Do NOT use external knowledge about the gene or mutation

QUESTION TYPES TO CONSIDER:
- Pathogenicity: "Is [mutation] pathogenic?"
- Disease association: "What disease/condition is associated with [mutation]?"
- Functional effect: "What is the functional effect of [mutation]?"
- Mechanism: "How does [mutation] affect protein function?"
- Clinical features: "What clinical features are seen with [mutation]?"
- Inheritance: "What is the inheritance pattern of [mutation]?"
- Frequency/prevalence: "How common is [mutation]?"

OUTPUT FORMAT (JSON only):
{
  "qa_pairs": [
    {
      "question": "...",
      "answer": "...",
      "rationale": "...",
      "evidence_sentence_ids": [list of ints],
      "question_type": "pathogenicity|disease|function|mechanism|clinical|inheritance|frequency|other"
    }
  ]
}

USER:
Generate 2-3 diverse QA pairs for this mutation:

MUTATION: {mutation}

EVIDENCE SENTENCES FROM THE PAPER:
{evidence_text}

Remember:
- Base every answer ONLY on the evidence sentences above
- Include "{mutation}" in every question AND answer
- Quote specific phrases from evidence in rationales
- Cite sentence numbers
- Explain how evidence supports the answer
- Try to cover different question types
- Include negative framing if appropriate

Output JSON only.
```

---

## 2. Independent answering (Phase 2)

`claim_extractor/phase_2/phase2_generator.py` · `deepseek/deepseek-chat`

**Inputs:** mutation, the Phase 1 *question*, and article context. This model
never sees Phase 1's answer or rationale — that isolation is what makes the
subsequent cross-validation meaningful.
**Output:** `answer`, `rationale`, `evidence_sentence_ids`, `confidence`,
`answerable`.

```
SYSTEM:
You are a scientific QA system. Answer the question about the specified mutation
using ONLY the provided context.

RULES:
1. Use ONLY information from the context - no external knowledge
2. If the answer is in the context, provide it - even if stated indirectly or as
   part of a group
3. Quote the relevant text and cite sentence numbers [N]
4. Standard mutation notation like "c.1033C>T, p.(Arg345Trp)" means c.1033C>T causes
   the p.Arg345Trp amino acid change

OUTPUT FORMAT (JSON only):
{
  "answer": "Your answer including the mutation name",
  "rationale": "Quote from context [N]. Brief explanation.",
  "evidence_sentence_ids": [N],
  "confidence": "high|medium|low",
  "answerable": true
}

Only set "answerable": false if the context truly contains NO relevant information
about the mutation.

USER:
Answer this question about the mutation {mutation}:

QUESTION: {question}

CONTEXT FROM RESEARCH PAPER:
{context}

Remember:
- Use ONLY the information in the context above
- Quote specific phrases and cite sentence numbers [N]
- Include "{mutation}" in your answer
- If the context doesn't contain enough information, say so

Output JSON only.
```

---

## 3. Cross-validation judge (Phase 3)

`claim_extractor/phase_3/phase3_judge.py` · `meta-llama/llama-3.1-8b-instruct`

**Inputs:** mutation, question, answer, rationale. **Output:** binary
`supported`. Run in both directions — Phase 1's rationale against Phase 2's
answer (`r1_supports_a2`) and vice versa (`r2_supports_a1`). A pair is labelled
`cross-grounded` only if both directions pass.

```
SYSTEM:
You are evaluating whether a rationale supports an answer to a question about a
genetic mutation.

Your task: Determine if the provided RATIONALE logically supports the provided
ANSWER to the QUESTION.

IMPORTANT:
- Focus on whether the rationale's evidence supports the answer's conclusion
- Check that Yes/No answers match the rationale's implications
- A rationale that describes the opposite of what the answer claims means NOT SUPPORTED
- The rationale must provide evidence for the specific claim in the answer

Output JSON only:
{
  "supported": true,
  "reasoning": "Brief explanation"
}

or

{
  "supported": false,
  "reasoning": "Brief explanation"
}

USER:
MUTATION: {mutation}

QUESTION: {question}

ANSWER: {answer}

RATIONALE: {rationale}

Does this rationale support this answer? Output JSON only.
```

---

## 4. QA-pair auditing judge (dataset quality)

`human_annotation/llm_annotation/run_llm_judges.py` ·
`google/gemini-2.0-flash-001` + `anthropic/claude-3.5-haiku`

Used to audit dataset quality against human annotators. **Inputs:** mutation,
question, answer, rationale. **Output:** 0–2 correctness + brief reasoning.

```
SYSTEM:
You are an expert scientific reviewer evaluating the correctness of answers about
genetic mutations.

You will be given:
1. A question about a genetic mutation
2. An answer to evaluate
3. A rationale/evidence supporting the answer

Rate the CORRECTNESS of the answer on this scale:
- 0 = Incorrect: The answer is wrong, contradicts the evidence, or makes unsupported claims
- 1 = Partially correct: The answer has some correct elements but is incomplete,
      imprecise, or has minor errors
- 2 = Completely correct: The answer is accurate, complete, and fully supported by
      the rationale

Output ONLY a JSON object with this format:
{
  "correctness": <0, 1, or 2>,
  "reasoning": "<brief explanation for your rating>"
}

USER:
Evaluate this answer about the mutation: {mutation}

QUESTION: {question}

ANSWER TO EVALUATE:
{answer}

RATIONALE/EVIDENCE:
{rationale}

Rate the correctness (0=incorrect, 1=partial, 2=correct). Output JSON only.
```

---

## 5. MutaPLM evaluation judge (UC1 / UC2)

`macaulay_eval/scripts/llm_judge_multi.py` · `deepseek/deepseek-v3.2` +
`google/gemini-2.0-flash-001` + `openai/gpt-4o-mini` · `temperature = 0`,
`max_tokens = 200`

Produces the "Judge mean" row. **Inputs:** exactly three fields — question,
ground-truth answer, prediction. The judge sees no gene, sequence, model identity
or system label. **Output:** three integer 0–5 scores. The reported figure is the
mean over the three dimensions, then over the three judges.

- **correctness** — factual accuracy of the predicted information
- **completeness** — how much of the ground-truth content is covered
- **relevance** — how directly the prediction addresses the question

```
You are an expert evaluator for protein-mutation Q&A.
Evaluate the predicted answer against the ground-truth answer.

Question: {question}

Ground Truth: {ground_truth}

Prediction: {prediction}

Score each on a 0-5 integer scale:
- correctness: factual accuracy of the predicted information
- completeness: how much of the ground-truth content is covered
- relevance: how directly the prediction addresses the question

Respond with ONLY a JSON object, no extra text:
{"correctness": 0-5, "completeness": 0-5, "relevance": 0-5}
```

### What fills `{question}` in each use case

The template above is shared by both use cases, but the question slot is filled
differently, and this changes what `relevance` measures. The two use cases are
specified as:

| UC | Input | Output | Judged for |
|:--|:--|:--|:--|
| UC1 | WT sequence + mutation | effect description | adequacy |
| UC2 | WT sequence + mutation + **question** | answer to the question | correctness / **relevance to that question** |

and the `{question}` slot is filled differently in each:

| | `{question}` is filled with | Consequence |
|:--|:--|:--|
| **UC1** | a per-item instruction naming the mutation (200/200): `Describe the effect of the A270D mutation in this protein.` | Correct as specified — effect description *is* the UC1 task, so the judge scores UC1 adequacy. |
| **UC2** | the constant **`Describe the effect of this mutation.`** — UC2 prediction files carry no `question` field (0/200 rows) | **The UC2 judge never saw the question and scored the answers under the UC1 criterion.** `relevance` measured relevance to effect description, not to the question. |

So the prompt every UC2 item was actually judged against is:

```
You are an expert evaluator for protein-mutation Q&A.
Evaluate the predicted answer against the ground-truth answer.

Question: Describe the effect of this mutation.

Ground Truth: {ground_truth}

Prediction: {prediction}

Score each on a 0-5 integer scale:
- correctness: factual accuracy of the predicted information
- completeness: how much of the ground-truth content is covered
- relevance: how directly the prediction addresses the question

Respond with ONLY a JSON object, no extra text:
{"correctness": 0-5, "completeness": 0-5, "relevance": 0-5}
```

`correctness` and `completeness` are scored against `{ground_truth}` — the item's
real reference answer — in both use cases, and are unaffected. Only `relevance`
is affected. UC2 correctness alone gives +0.92 to +1.35 across the three judges,
against +1.37 for the three-dimension mean, so the direction and approximate
magnitude of the reported UC2 result do not depend on the affected dimension.

Two related artefacts, for anyone reproducing this:

- `judge_inputs/qa_1500_*_nondelex.json` carries a `question` field built from
  `mutation_type`, yielding `Describe the effect of the Substitution mutation in
  this protein.` It is degenerate for the same reason and should not be used.
- The UC2 items' real questions are recoverable — all 200 join to
  `qa_test_1500.json` on the reference answer (199 uniquely).

---

## 6. Retrieval-experiment judge

`rag/judge.py` · `google/gemini-2.5-flash-lite` + `anthropic/claude-haiku-4.5` +
`deepseek/deepseek-v3.2` · `temperature = 0`, `seed = 42`

Grades the leave-one-out retrieval experiment on the same 0/1/2 rubric as role 4,
so the two are directly comparable. **Inputs:** question, reference answer,
prediction — the judge is blind to which arm produced the prediction. **Output:**
`score` plus a separate `abstained` flag, so a declined answer is not conflated
with a confident error.

No OpenAI judge is used here: the generator set includes `openai/gpt-5.6-terra`
and same-family self-preference is a known judge bias. The strict (`score == 2`)
and lenient (`score >= 1`) decision rules were fixed in code before any results
existed.

```
SYSTEM:
You are grading answers about genetic variants against a reference answer. You grade
factual agreement only.

USER:
Question: {question}

Reference answer: {gold}

Predicted answer: {pred}

Grade the predicted answer against the reference on this scale:
  2 = fully correct - asserts the same fact(s) as the reference, no contradictions
  1 = partially correct - some reference content is right, but it is incomplete,
      vague, or mixes correct and incorrect claims
  0 = incorrect - contradicts the reference, or asserts something the reference
      does not support

Grade FACTUAL AGREEMENT ONLY. Ignore length, phrasing, fluency, confidence, and
any stylistic resemblance to the reference wording. A blunt one-clause answer
that agrees factually scores 2; a fluent, well-written answer that disagrees
scores 0. Do not reward wording that merely sounds like the reference.

If the prediction declines to answer (for example "INSUFFICIENT INFORMATION"),
set abstained=true and score=0.

Respond with ONLY this JSON object and nothing else:
{"score": 0, "abstained": false}
```
