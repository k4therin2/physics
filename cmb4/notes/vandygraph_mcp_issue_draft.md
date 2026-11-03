# Vandygraph MCP server: gaps hit while publishing a page from Claude Code (draft issue)

Context: publishing and updating an ~80 KB single-file HTML report from a Claude Code session via the Vandygraph MCP tools
(create_upload, replace_upload_html, append_upload_html, create_site_upload, set_upload_file, delete_upload, get_upload, rename_upload).

## 1. No way to change visibility of an existing upload
`create_upload` accepts `visibility`, but there is no `update_upload` / `set_visibility` tool and `replace_upload_html` does not take
visibility. To make a private page public I had to `delete_upload` and `create_upload` again (the slug happened to be reusable, so
the share URL survived, but that is luck). Request: `set_upload_visibility(id, visibility)` or a `visibility` field on
`replace_upload_html` / `rename_upload`.

## 2. Content must be passed inline; no file path or URL source
Every tool takes the full document as an inline string. For an agent, an 80 KB page is ~40k tokens of output per publish, and
the Read tool caps at ~25k tokens, so the file had to be split and sent as `replace_upload_html` + `append_upload_html`.
Request: accept a local file path (for MCP servers running with filesystem access) or a fetchable URL / multipart upload, and/or
return a presigned PUT URL so the agent can `curl` the bytes without routing them through the model.

## 3. No binary assets for site uploads
`create_site_upload` / `set_upload_file` are text-only, so screenshots or images cannot be included in a multi-page site.
Request: binary file support (base64 or presigned upload).

## 4. Minor
- `get_upload` metadata does not include byte size / last-modified in the create response (only `byte_size` on replace).
- `list_uploads` would benefit from `workspace` + `visibility` filters.
- A `preview` / render-check endpoint (HEAD of the published page returning status and size) would let agents verify a publish.

Environment: Claude Code (claude.ai cloud session), Vandygraph MCP connector, 2026-11-03.
