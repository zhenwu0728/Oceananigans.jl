module Logger

export ModelLogger, Diagnostic, Setup

using Dates
using Logging

#####
##### Custom LogLevels
#####

_custom_log_level_docs = """
    Severity Order:
    Debug < Diagnostic < Setup < Info

    Usage:
    @logmsg Logging.LogLevel "Log Message"

    @logmsg comes from Base/Logging
    LogLevel can be any Base/Logging.LogLevel
    Log Message can be any expression that evaluates to a string (preferably human readable!)
"""

const Diagnostic = Logging.LogLevel(-500)  # Sits between Debug and Info
const Setup      = Logging.LogLevel(-125)

#####
##### ModelLogger
#####

"""
    ModelLogger(stream::IO, level::LogLevel)

Based on Logging.SimpleLogger it tries to log all messages in the following format

    message --- [dd/mm/yyyy HH:MM:SS] log_level source_file:line_number

The logger will handle any message from Diagnostic up by default.
"""
struct ModelLogger <: Logging.AbstractLogger
            stream :: IO
         min_level :: Logging.LogLevel
    message_limits :: Dict{Any,Int}
end

ModelLogger(stream::IO=stderr, level=Diagnostic) = ModelLogger(stream, level, Dict{Any,Int}())

Logging.shouldlog(logger::ModelLogger, level, _module, group, id) = get(logger.message_limits, id, 1) > 0

Logging.min_enabled_level(logger::ModelLogger) = logger.min_level

Logging.catch_exceptions(logger::ModelLogger) = false

function level_to_string(level::Logging.LogLevel)
    level == Diagnostic   && return "Diagnostic"
    level == Setup        && return "Setup"
    level == Logging.Warn && return "Warning"
    return string(level)
end

function Logging.handle_message(logger::ModelLogger, level, message, _module, group, id, filepath, line; maxlog = nothing, kwargs...)
    if maxlog !== nothing && maxlog isa Integer
        remaining = get!(logger.message_limits, id, maxlog)
        logger.message_limits[id] = remaining - 1
        remaining > 0 || return
    end

    buf = IOBuffer()
    iob = IOContext(buf, logger.stream)
    level_name = level_to_string(level)

    module_name = something(_module, "nothing")
    file_name   = something(filepath, "nothing")
    line_number = something(line, "nothing")
    msg_timestamp = Dates.format(Dates.now(), "[dd/mm/yyyy HH:MM:SS]")

    formatted_message = "$msg_timestamp $message ---  $level_name $file_name:$line_number"

    println(iob, formatted_message)
    write(logger.stream, take!(buf))

    return nothing
end

end
