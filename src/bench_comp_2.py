#!/usr/bin/env python
import csv
import sys
from typing import List, NamedTuple
from collections import defaultdict, Counter

class ToolData(NamedTuple):
    command: str
    mean: float
    stddev: float
    median: float
    user: float
    system: float
    min: float
    max: float

def parse_csv(file_path: str) -> List[ToolData]:
    data = []
    with open(file_path, 'r') as csvfile:
        try:
            reader = csv.reader(csvfile)
        except csv.Error as e:
            sys.exit('file {}, line {}: {}'.format(file_path, reader.line_num, e))
        
        try:
            next(reader)  # Skip header
        except StopIteration:
            return data
        for row in reader:
            command, mean, stddev, median, user, system, min_time, max_time = row
            tool = command.split()[0]
            # if tool is a path, extract the basename
            if '/' in tool:
                tool = tool.split('/')[-1]
            data.append(ToolData(
                tool,
                float(mean),
                float(stddev),
                float(median),
                float(user),
                float(system),
                float(min_time),
                float(max_time)
            ))
    return data

def process_files(file_paths: List[str]) -> List[ToolData]:
    all_data = []
    for file_path in file_paths:
        all_data.extend(parse_csv(file_path))
    return all_data

def rank_tools(data: List[ToolData]) -> List[ToolData]:
    # Group by command and calculate average mean for each command
    command_data = defaultdict(list)
    for tool in data:
        command_data[tool.command].append(tool)
    
    # Calculate average ToolData for each command
    averaged_data = []
    for command, tools in command_data.items():
        avg_tool = ToolData(
            command=command,
            mean=sum(t.mean for t in tools) / len(tools),
            stddev=sum(t.stddev for t in tools) / len(tools),
            median=sum(t.median for t in tools) / len(tools),
            user=sum(t.user for t in tools) / len(tools),
            system=sum(t.system for t in tools) / len(tools),
            min=min(t.min for t in tools),
            max=max(t.max for t in tools)
        )
        averaged_data.append(avg_tool)
    
    # Sort by mean execution time
    return sorted(averaged_data, key=lambda x: x.mean)

def count_fastest_tools(file_paths: List[str]) -> Counter:
    fastest_tools = Counter()
    for file_path in file_paths:
        data = parse_csv(file_path)
        if data:
            fastest_tool = min(data, key=lambda x: x.mean).command
            fastest_tools[fastest_tool] += 1
    return fastest_tools

def print_ranked_tools(ranked_tools: List[ToolData]):
    print("Ranked list of tools (from fastest to slowest):")
    print("-" * 80)
    print(f"{'Rank':<5}{'Command':<20}{'Mean (s)':<12}{'Median (s)':<12}{'StdDev (s)':<12}{'Min (s)':<12}{'Max (s)':<12}")
    print("-" * 80)
    for i, tool in enumerate(ranked_tools, 1):
        print(f"{i:<5}{tool.command:<20}{tool.mean:<12.6f}{tool.median:<12.6f}{tool.stddev:<12.6f}{tool.min:<12.6f}{tool.max:<12.6f}")

def print_fastest_tool_frequency(fastest_tools: Counter):
    print("\nRanked list of tools by frequency of being fastest:")
    print("-" * 50)
    print(f"{'Rank':<5}{'Command':<20}{'Frequency':<10}")
    print("-" * 50)
    for i, (tool, count) in enumerate(fastest_tools.most_common(), 1):
        print(f"{i:<5}{tool:<20}{count:<10}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <csv_file1> <csv_file2> ...")
        sys.exit(1)
    
    file_paths = sys.argv[1:]
    all_data = process_files(file_paths)
    ranked_tools = rank_tools(all_data)
    
    print_ranked_tools(ranked_tools)
    
    fastest_tools = count_fastest_tools(file_paths)
    print_fastest_tool_frequency(fastest_tools)

if __name__ == "__main__":
    main()