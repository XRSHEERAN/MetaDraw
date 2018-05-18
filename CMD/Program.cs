using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Text.RegularExpressions;

namespace CMD
{
    class Program
    {
        static void Main(string[] args)
        {
            //setup
            if (false)
            {
                List<string> commands = new List<string>();
                commands.Add("echo \"Checking for updates and installing any missing dependencies. Please enter your password for this step:\n\"");
                commands.Add("sudo easy_install pip");
                commands.Add("sudo -H pip install --upgrade pip numpy scipy");

                string installationLogsDirectory = Path.Combine((Assembly.GetEntryAssembly().Location), "scripts", "installLogs");
                Directory.CreateDirectory(installationLogsDirectory);
            }

            //Run pdeep python script 
            List<string> command = new List<string>();
            string workingPath = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location));
            string scriptPath = Path.Combine(workingPath, "test.bash");
            command.Add("cd" + " " + ConvertWindowsPath(workingPath));
            command.Add("python3 example.py > test.txt");
            GenerateAndRunScript(scriptPath, command).WaitForExit();
        }

        public static Process GenerateAndRunScript(string scriptPath, List<string> commands)
        {
            GenerateScript(scriptPath, commands);
            return RunBashCommand(@"bash", ConvertWindowsPath(scriptPath));
        }

        public static void GenerateScript(string scriptPath, List<string> commands)
        {
            Directory.CreateDirectory(Path.GetDirectoryName(scriptPath));
            using (StreamWriter writer = new StreamWriter(scriptPath))
            {
                foreach (string cmd in commands)
                {
                    writer.Write(cmd + "\n");
                }
            }
        }

        public static Process RunBashCommand(string command, string arguments)
        {
            Process proc = new Process();
            proc.StartInfo.FileName = @"C:\Windows\System32\bash.exe";
            proc.StartInfo.Arguments = "-c \"" + command + " " + arguments + "\"";
            proc.Start();
            return proc;
        }

        private static Regex driveName = new Regex(@"([A-Z]:)");
        private static Regex forwardSlashes = new Regex(@"(\\)");

        public static string ConvertWindowsPath(string path)
        {
            if (path == null) return null;
            if (path == "") return "";
            if (path.StartsWith("/mnt/")) return path;
            return "/mnt/" + Char.ToLowerInvariant(path[0]) + driveName.Replace(forwardSlashes.Replace(path, "/"), "");
        }

    }
}
