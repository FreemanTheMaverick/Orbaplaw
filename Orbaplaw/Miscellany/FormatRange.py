import itertools

def FormatRange(numbers):
	if not numbers:
		return ""
	sorted_unique = sorted(set(numbers))
	groups = []
	for key, group in itertools.groupby(enumerate(sorted_unique), lambda x: x[1] - x[0]):
		group = list(group)
		start = group[0][1]
		end = group[-1][1]
		groups.append((start, end))
	parts = []
	for s, e in groups:
		parts.append(f"{s}-{e}" if s != e else f"{s}")
	return ",".join(parts)
