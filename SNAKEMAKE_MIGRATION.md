# Snakemake 7 → 9 Migration Guide

## Version Inconsistency Resolution

### Current State
- **README.md:** Claims Snakemake 5.22 (outdated)
- **nccrPipe_environment.yaml:** Specifies `snakemake>=7.1.*`
- **Actual Version:** Likely 7.x (need to verify with users)

### Target State
- **Python:** 3.12.*
- **Snakemake:** 9.x (latest stable)

---

## Breaking Changes: Snakemake 7 → 9

### 1. Conda Environment Handling

#### Old (Snakemake 7)
```python
rule example:
    conda: "envs/myenv.yaml"  # Can use relative path
```

#### New (Snakemake 9)
```python
rule example:
    conda: "envs/myenv.yaml"  # Must be relative to Snakefile or use workflow.source_path
```

**Action Required:**
- ✅ Current code uses relative paths - should work
- ⚠️ Need to test all conda environments
- Consider using `workflow.source_path()` for clarity

---

### 2. Output Markers with `touch()`

#### Old (Snakemake 7)
```python
output: marker = touch("file.done")
```

#### New (Snakemake 9)
```python
# touch() is deprecated
output: marker = "file.done"
shell: "touch {output.marker}"
```

**Action Required:**
- 🚨 **HIGH PRIORITY** - `touch()` used extensively throughout codebase
- Examples in `call_variants.smk:21`, `call_variants.smk:49`, etc.
- Need to replace all instances

**Files affected:**
- `workflow/rules/call_variants.smk`
- `workflow/rules/assemble.smk`
- `workflow/rules/annotate.smk`
- `workflow/rules/preprocess.smk`
- And others...

---

### 3. Cluster Execution

#### Old (Snakemake 7)
```bash
snakemake --cluster "sbatch {params}" --jobs 10
```

#### New (Snakemake 9)
```bash
# Profiles are recommended
snakemake --profile profiles/slurm --jobs 10

# Or use executor
snakemake --executor slurm --jobs 10
```

**Action Required:**
- 🔄 **MAJOR REFACTOR** - Current implementation in `main.py:233-256`
- Hardcoded cluster submission strings
- Should migrate to Snakemake profiles
- Create `profiles/slurm/config.yaml`

---

### 4. Resource Specification

#### Old (Snakemake 7)
```python
params:
    mem=7700,
    time=1400,
    qoutfile="log.out",
    qerrfile="log.err"
```

#### New (Snakemake 9)
```python
# Use resources instead of params for cluster resources
resources:
    mem_mb=7700,
    runtime=1400

# Log files handled by profile
log: "logs/rule.log"
```

**Action Required:**
- 🚨 **MAJOR CHANGE** - All rules use `params:` for resources
- Need to migrate to `resources:` directive
- Update cluster profile to handle resources

**Pattern to replace:**
```python
# FROM:
params:
    qerrfile='path.qerr',
    qoutfile='path.qout',
    scratch=6000,
    mem=7700,
    time=1400

# TO:
resources:
    mem_mb=7700,
    runtime=1400,
    disk_mb=6000
log:
    'path.log'
```

---

### 5. Config File Handling

#### Changes
- `yaml.load()` deprecated → use `yaml.safe_load()`
- Config validation recommended

**Action Required:**
- 🚨 **SECURITY ISSUE** - `workflow/main.py:42` uses `yaml.load()`
- Replace with `yaml.safe_load()`
- Add config schema validation

---

### 6. Shell Command Syntax

#### Old (Snakemake 7)
```python
shell: "command {input} > {output}"
```

#### New (Snakemake 9)
```python
# More strict about wildcards and escaping
shell: "command {input} > {output}"
# Use f-strings in Python code, not in shell:
```

**Action Required:**
- ⚠️ Review all shell commands for syntax issues
- Test with Snakemake 9

---

### 7. Python Version Requirements

#### Changes
- Python 3.7 deprecated
- Python 3.8+ required for Snakemake 8
- Python 3.9+ recommended for Snakemake 9
- Python 3.12 supported

**Action Required:**
- ✅ Update to Python 3.12
- Test all Python scripts for 3.12 compatibility

---

## Migration Strategy

### Phase 1: Preparation
1. ✅ Create refactor branch
2. ✅ Document current state
3. ✅ Set up testing
4. 🔄 Audit all uses of deprecated features

### Phase 2: Critical Fixes (Do First)
1. **Replace `touch()` outputs**
   - Search: `touch\(`
   - Replace with explicit shell command
   - Priority: HIGH (breaks in Snakemake 9)

2. **Fix `yaml.load()` security issue**
   - File: `workflow/main.py:42`
   - Replace with `yaml.safe_load()`
   - Priority: HIGH (security)

3. **Update resource specifications**
   - Migrate `params:` → `resources:`
   - All rule files affected
   - Priority: HIGH (required for profiles)

### Phase 3: Cluster Modernization
1. **Create Snakemake profiles**
   - `profiles/slurm/config.yaml`
   - `profiles/local/config.yaml`
   - Move cluster logic from Python to profiles

2. **Update CLI**
   - Remove cluster string building
   - Use `--profile` flag
   - Update documentation

### Phase 4: Testing
1. Test with Snakemake 7.x (current)
2. Test with Snakemake 8.x (intermediate)
3. Test with Snakemake 9.x (target)

---

## Detailed Audit Results

### `touch()` Usage Count
```bash
# Run this to find all instances:
grep -r "touch(" workflow/rules/*.smk
```

**Found in:**
- `call_variants.smk`: ~10 instances
- `assemble.smk`: ~5 instances
- `annotate.smk`: ~3 instances
- `preprocess.smk`: ~4 instances
- `compare_genomes.smk`: ~2 instances
- Others: TBD

### Cluster Resource Params Count
**All rules have:**
```python
params:
    qerrfile=...,
    qoutfile=...,
    scratch=...,
    mem=...,
    time=...
```

**Estimated:** ~50+ rules to update

---

## Backward Compatibility Notes

### For Existing Users

**Option 1: Pin to Snakemake 7**
```yaml
# environment.yaml
snakemake=7.*
python=3.8.*
```
- Keep using current codebase
- No changes needed
- Stick with master branch

**Option 2: Upgrade to Snakemake 9**
```yaml
# environment.yaml
snakemake>=9.0
python>=3.12
```
- Use refactored codebase
- Follow migration guide
- Switch to refactor branch

### Migration Path for Users

1. **Current users:** Stay on master with Snakemake 7
2. **New users:** Use refactor branch with Snakemake 9
3. **Transition period:** Both versions maintained
4. **After release:** v1.0 becomes default

---

## Testing Strategy

### Test Matrix
| Python | Snakemake | Status |
|--------|-----------|--------|
| 3.8    | 7.x       | Current working |
| 3.12   | 7.x       | Test compatibility |
| 3.12   | 8.x       | Intermediate test |
| 3.12   | 9.x       | Target |

### Test Workflows
1. Variant calling (most used)
2. Assembly
3. Preprocessing only
4. Full pipeline

### Test Modes
1. Dry run (`-np`)
2. Local execution
3. Cluster execution (if available)

---

## Timeline

### Week 1-2: Critical Fixes
- Remove `touch()`
- Fix `yaml.load()`
- Create Snakemake profiles

### Week 3-4: Resource Migration
- Update all rules
- Migrate params to resources
- Test incrementally

### Week 5-6: Testing
- Test matrix execution
- Fix issues
- Documentation

---

## Resources

### Snakemake Documentation
- [Snakemake 9 Changelog](https://snakemake.readthedocs.io/en/stable/project_info/history.html)
- [Cluster Execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
- [Profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)

### Migration Examples
- [Snakemake profile examples](https://github.com/Snakemake-Profiles)
- [SLURM profile](https://github.com/Snakemake-Profiles/slurm)

---

## Next Steps

1. ✅ Document complete
2. 🔜 Audit all `touch()` usage
3. 🔜 Create SLURM profile template
4. 🔜 Test one rule file migration
5. 🔜 Scale to all rules
