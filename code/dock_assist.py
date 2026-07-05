import os
import glob
from PIL import Image
from collections import Counter
from typing import Annotated, TypedDict
import time, sys
from ollama import Client as ollama_client
from rich.console import Console
from rich.markdown import Markdown

from modrag_protein_functions import get_protein_from_pdb, find_PDBID_node, smiles_node
from vina_dock import blind_dock_agent
console = Console(width=80)
import modrag_protein_functions


# get the print flag from the command line arguments, if they exist, otherwise set it to False
if len(sys.argv) > 1 and sys.argv[1] == '--print':
    print_flag = True
else:
    print_flag = False

# Set print_flag in all imported modules
modrag_protein_functions.print_flag = print_flag

tools = [find_PDBID_node, get_protein_from_pdb, smiles_node, blind_dock_agent]

#get ket from shell variable $OLLAMA_API_KEY
ollama_key = os.getenv('OLLAMA_API_KEY')

# models = ['deepseek-v3.1:671b', 'gpt-oss:120b', 'gpt-oss:20b',
#           'devstral-2:123b', 'cogito-2.1:671b',
#           'nemotron-3-nano:30b', 'gemini-3-flash-preview',
#           'kimi-k2:1t', 'kimi-k2.5', 'gemma4:31b-cloud']

models = [
    'gemma4:31b', 'glm-5.2', 'kimi-k2.7-code',
    'deepseek-v4-pro', 'qwen3.5:397b',
]

model = models[0] # default model to use for chat

sys_message = f'''
You are a drug discovery assistant named Modrag-dock. You have access to the 
following tools: {', '.join([tool.__name__ for tool in tools])}. Your job is 
to find and retreive pdb files for proteins of interest and then perform docking simulations. 
You will be given a protein name and you will search the protein databank for PDB IDs that match 
along with the entry titles. You will pick the most relevant PDB ID and use 
it to retreive the pdb file for that protein. You will then dock any requested ligands to the 
protein and return the results.
'''

global messages
messages = [{'role': 'system', 'content': sys_message}]

def start_chat():
  '''
  Initializes a new chat session by resetting the chat history, reasoning, and messages.
  '''
  global chat_history, messages, reasoning
  chat_history = []
  reasoning = []
  messages = [{'role': 'system', 'content': sys_message}]

def chat_turn(prompt: str):
  '''
  Handles a single turn of the chat by sending the user's prompt to the Ollama API,
  processing the response, and executing any tool calls if present.
  '''
  global chat_history, messages, reasoning
  
  client = ollama_client(host = 'https://ollama.com',
            headers={'Authorization': f'Bearer {ollama_key}'})

  available_functions = {
    'find_PDBID_node': find_PDBID_node,
    'get_protein_from_pdb': get_protein_from_pdb,
    'smiles_node': smiles_node,
    'blind_dock_agent': blind_dock_agent
  }

  messages.append({'role': 'user', 'content': prompt})

  while True:
      response = client.chat(
          model=model,
          messages=messages,
          tools=[find_PDBID_node, get_protein_from_pdb, smiles_node, blind_dock_agent],
          think=True,
      )
      messages.append(response.message)
      if print_flag:
          print('------------------------------------------------------------------------')
          print("Thinking: ", response.message.thinking)
          print('------------------------------------------------------------------------')
          print("Content: ", response.message.content)
          print('------------------------------------------------------------------------')
      if response.message.tool_calls:
        for tc in response.message.tool_calls:
          if tc.function.name in available_functions:
            if print_flag:
              print(f"Calling {tc.function.name} with arguments {tc.function.arguments}")
            result = available_functions[tc.function.name](**tc.function.arguments)
            if print_flag:
              print(f"Result: {result}")
              print('------------------------------------------------------------------------')
            # add the tool result to the messages
            messages.append({'role': 'tool', 'tool_name': tc.function.name, 'content': str(result)})
      else:
        # end the loop when there are no more tool calls
        break

  return '', None, messages[-1]['content']

start_chat()

# Clean up any .sdf files left in pdb_files/ from a previous session so the
# next docking run starts fresh. .pdb files (the prepared receptors) are kept.
sdf_files = glob.glob('pdb_files/*.sdf')
for sdf in sdf_files:
  try:
    os.remove(sdf)
  except OSError as e:
    print(f"\033[38;5;208mWarning: could not remove {sdf}: {e}\033[0m")
if sdf_files:
  print(f"\033[38;5;208mNote: removed {len(sdf_files)} leftover .sdf file(s) from pdb_files/.\033[0m")

header_string = f'''
\033[1;36m*******************************************************
*                                                     *
\033[1;35m*  __  __       ____        _                         *
\033[1;36m* |  \/  | ___ |  _ \ _ __ / \   __ _                 *
\033[1;35m* | |\/| |/ _ \| | | | '__/ _ \ / _` |                *
\033[1;36m* | |  | | (_) | |_| | | / ___ \ (_| |                *
\033[1;35m* |_|__|_|\___/|____/|_|/_/   \_\__, |     _     _    *
\033[1;36m* |  _ \  ___   ___| | __    / \|___/_ ___(_)___| |_  *
\033[1;35m* | | | |/ _ \ / __| |/ /   / _ \ / __/ __| / __| __| *
\033[1;36m* | |_| | (_) | (__|   <   / ___ \\__ \__ \ \__ \ |_  *
\033[1;35m* |____/ \___/ \___|_|\_\ /_/   \_\___/___/_|___/\__| *
\033[1;36m*                                                     *
*******************************************************\033[0m
\033[1;31mThe MOdular DRug design AGent!\033[0m
\033[1;35mA command-line interface (CLI) for drug\033[0m
\033[1;31mdiscovery and molecular design.\033[0m
\033[0m'''


print(header_string)

next_prompt = input("\033[1;36mI'm the MoDrAg docking assistant! \
I can help with protein docking.\nI just need a protein name and \
names or smiles for the ligands.\nWhat can I help with today? > \033[0m")
print('')
if next_prompt == 'quit':
  print("\033[1;35mResponse > \033[0mGoodbye!")
else:
  start_time = time.time()
  _, _, response_content = chat_turn(next_prompt)
  end_time = time.time()

  time_for_inf = (end_time - start_time) / 60
  print(f"\033[1;35mResponse {time_for_inf:.2f}m > \033[0m")
  console.print(Markdown(response_content))

  print(f"\033[38;5;208mNote: The .pdb file used in the docking and the docked .sdf pose files "
        f"can be found in the pdb_files/ directory. You can view them in PyMOL.\033[0m")

while next_prompt != 'quit':
  print('')
  next_prompt = input("\033[1;36mWhat else can I help with? > \033[0m")
  print('')
  if next_prompt == 'quit':
    print("\033[1;35mResponse > \033[0mBe sure to remove any SDF files you need before next time, \n\
as they will be cleaned at the next start-up. Goodbye!")
    break
  else:
    start_time = time.time()
    _, _, response_content = chat_turn(next_prompt)
    end_time = time.time()

    time_for_inf = (end_time - start_time) / 60
    print('')
    print(f"\033[1;35mResponse {time_for_inf:.2f}m > \033[0m")
    console.print(Markdown(response_content))

    print(f"\033[38;5;208mNote: The .pdb file used in the docking and the docked .sdf pose files "
          f"can be found in the pdb_files/ directory. You can view them in PyMOL.\033[0m")